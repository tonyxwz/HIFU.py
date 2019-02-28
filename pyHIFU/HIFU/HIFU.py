from multiprocessing import current_process, Pool
import numpy as np

from pyHIFU.physics.medium import InitMedium, Medium, MediaComplex
from pyHIFU.transducer import TElement, Transducer
from pyHIFU.box import Sparse3D, Box


class HIFU(object):
    """ HIFU system (methods encapsulated)
    """

    def __init__(self, transducer: Transducer):
        self.transducer = transducer

    def set_transducer(self,
                       nature_focus=0,
                       actual_focus=0,
                       focus_diameter=0,
                       frequency=0,
                       element_coordinates=None,
                       element_radius=None,
                       element_power=None,
                       element_properties=None,
                       **kw):
        self.transducer = Transducer(nature_focus, actual_focus,
                                     focus_diameter, frequency,
                                     element_coordinates, element_radius,
                                     element_power, element_properties, **kw)

    def run(self,
            init_medium: InitMedium,
            mc: MediaComplex,
            box: Box,
            n_rays=500,
            trident_angle=5e-3,
            theta_max=0.1,
            n_core=4,
            verbose=True):
        complex_pressure = np.zeros(box.nxyz, dtype=np.complex128)
        if n_core > 1:
            pool = Pool(n_core)
            async_results = []
            for i, co in enumerate(self.transducer.element_coordinates):
                ar = pool.apply_async(
                    self.run_te_wrapper,
                    args=(
                        i,
                        co,
                        self.transducer.element_radius,
                        self.transducer.element_power,
                        self.transducer.frequency,
                        self.transducer.nature_focus,
                        init_medium,
                        0,  # init phase
                        n_rays,
                        trident_angle,
                        theta_max,
                        mc,  # media complex
                        box,  # sampling box
                        verbose,  # print message to stderr
                        True  # whether running in pararllel
                    ))
                async_results.append(ar)

            pool.close()
            pool.join()

            for ar in async_results:
                complex_pressure += ar.get()
        else:
            for i, co in enumerate(self.transducer.element_coordinates):
                complex_pressure += self.run_te_wrapper(
                    i,
                    co,
                    self.transducer.element_radius,
                    self.transducer.element_power,
                    self.transducer.frequency,
                    self.transducer.nature_focus,
                    init_medium,
                    0,  # init phase
                    n_rays,
                    trident_angle,
                    theta_max,
                    mc,  # media complex
                    box,  # sampling box
                    verbose,  # print message to stderr
                    False  # whether running in pararllel
                )
        self.result = complex_pressure  # TODO to heat production
        return complex_pressure

    @staticmethod
    def run_te_wrapper(el_id,
                       center,
                       radius,
                       power,
                       freq,
                       nature_f,
                       init_medium,
                       initial_phase,
                       n_rays,
                       trident_angle,
                       theta_max,
                       mc,
                       box,
                       verbose,
                       parallel=False):
        if parallel:
            pname = current_process().name
        te = TElement(
            el_id,
            center,
            radius=radius,
            power=power,
            freq=freq,
            nature_f=nature_f)
        te.initialize(
            init_medium,
            initial_phase=initial_phase,
            n=n_rays,
            trident_angle=trident_angle,
            theta_max=theta_max)
        interface = init_medium.shape[0]
        for tr in te:
            tr.pow_ray.end = interface.intersect_line(tr.pow_ray)
            tr.aux_ray1.end = interface.intersect_line(tr.aux_ray1)
            tr.aux_ray2.end = interface.intersect_line(tr.aux_ray2)

        bundle_dict = te.cast(mc)

        pc = np.zeros(box.nxyz, dtype=np.complex128)

        for _bundle_str, tr_list in bundle_dict.items():
            I = Sparse3D(box.nxyz)
            ph = Sparse3D(box.nxyz)
            counter = Sparse3D(box.nxyz, dtype=int)
            for tr in tr_list:
                box.intersect_trident(tr, I, ph, counter)

            for k in I.getdata():
                # one could just assume I, ph, counter always have the same keys
                # use the Z of last tr because one bundle have the same medium
                pc[k] += np.sqrt(2 * tr.medium.Z * I[k] / counter[k]) * np.exp(
                    1j * ph[k] / counter[k])
        if verbose:
            if parallel:
                print(f"{pname}: TE# {el_id} processed")
            else:
                print(f"TE# {el_id} is processed")

        return pc
