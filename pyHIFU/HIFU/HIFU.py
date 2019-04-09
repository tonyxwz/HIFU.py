import time
from multiprocessing import Pool, current_process

import numpy as np

from pyHIFU.box import Box, Sparse3D
from pyHIFU.physics.medium import MediaComplex, Medium
from pyHIFU.transducer import TElement, Transducer


class HIFU(object):
    """ HIFU system (methods encapsulated)
    """

    def __init__(self, transducer: Transducer):
        self.transducer = transducer
        seed = int(time.time())
        seed = 34839194
        print("random seed:", seed)
        # np.random.seed(18973894)
        np.random.seed(34839194)
        self.result = None

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
            mc: MediaComplex,
            box: Box,
            nrays=500,
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
                        0,  # init phase
                        nrays,
                        trident_angle,
                        theta_max,
                        mc,  # media complex
                        box,  # sampling box
                        verbose,  # print message to stderr
                        True  # whether running in parallel
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
                    0,  # init phase
                    nrays,
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
                       initial_phase,
                       nrays,
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

        bundle_dict = te.cast(
            mc=mc,
            initial_phase=initial_phase,
            nrays=nrays,
            trident_angle=trident_angle,
            theta_max=theta_max)

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
                print(f"{pname}: TE #{el_id} processed")
            else:
                print(f"TE #{el_id} is processed")

        return pc
