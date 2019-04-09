import time
from collections import deque
from multiprocessing import Pool, current_process, Manager

import numpy as np
from cached_property import cached_property
from scipy import integrate, special

from .geometric.curves import Ring
from .geometric.lines import Ray
from .geometric.surfaces import Plane, Sphere
from .geometric.vec3 import Vec3
from .io.config import readjson
from .physics import LONGITUDINAL, SHEAR
from .ray import Trident


class TElement(list):
    """ transducer element class """

    def __init__(self,
                 el_id,
                 center,
                 radius=0,
                 power=0,
                 freq=0,
                 nature_f=np.array([0, 0, 0]),
                 **kw):
        """ only assgin necessary parameters here """
        super().__init__()  # len(self) == 0
        self.el_id = el_id
        self.center = center  # also the start of all tridents
        self.radius = radius
        self.required_power = power  # needed in ffa
        self.frequency = freq
        self.nature_f = nature_f
        self.axial_ray = Ray(self.center, self.nature_f - self.center)

    def fluxfunc(self, x):
        return np.sin(x) * (2 * special.jv(1, self.ka * x) / (self.ka * x))**2

    def initialize(self,
                   init_medium,
                   initial_phase=0,
                   n=100,
                   trident_angle=1e-4,
                   theta_max=np.pi / 6):
        """ use `te.cast` directly
        initialize all trident rays until they hit markoil interface
        `init_medium`: e.g. lossless / markoil
        `n`: number of rays per transducer
        `trident_angle`: the angle between powray and auxray
        `sampling_angle`: theta_max
        """
        self.trident_angle = trident_angle
        self.theta_max = theta_max
        AA = 1 - np.cos(theta_max)
        self.nrays = n
        self.init_medium = init_medium
        # self.fluxfunc = lambda x: np.sin(x) * (2 * special.jv(1, self.ka * x) / (self.ka * x) )**2
        self.initial_phase = initial_phase
        vr = self.axial_ray.perpendicularDirection()
        vr = Vec3.rotate(vr, self.axial_ray.d, np.random.random() * np.pi * 2)
        for i in range(self.nrays):
            # initialize n_rays number of random directed trident rays
            # theta = np.random.random() * self.theta_max
            theta = np.arccos(1 - AA * np.random.random())
            p1 = self.axial_ray.to_coordinate(np.cos(theta))

            # random angle on the ring by counter-clockwise rotation
            beta = np.random.random() * np.pi * 2
            v_ = Vec3.rotate(vr, self.axial_ray.d, beta)

            p_end = p1 + Vec3.normalize(v_) * np.sin(theta)
            pow_dire = p_end - self.center
            ray_helper = Ray(self.center, pow_dire)

            # two perpendicular directions for auxrays
            n_helper = ray_helper.perpendicularDirection()
            n2_helper = np.cross(pow_dire, n_helper)
            l = np.tan(self.trident_angle) * np.linalg.norm(pow_dire)
            assert l > 0
            a1_end = p_end + Vec3.normalize(n_helper) * l
            a2_end = p_end + Vec3.normalize(n2_helper) * l

            # calculate initial power at theta
            z = self.distance_z[LONGITUDINAL]
            pressure0 = self.ffa(z, theta)
            I0 = np.abs(pressure0)**2 / (2 * self.init_medium.Z[LONGITUDINAL])
            self.append(
                Trident(
                    self.center,
                    pow_dire,
                    self.center,
                    a1_end - self.center,
                    self.center,
                    a2_end - self.center,
                    I0,
                    self.frequency,
                    self.initial_phase,
                    len0=z,
                    el_id=self.el_id,
                    ray_id=self.el_id * self.nrays + i,
                    medium=init_medium,
                    legacy=[],
                    wave_type=LONGITUDINAL))

    def cast(self,
             nrays=100,
             mc=[],
             trident_angle=1e-4,
             theta_max=np.pi / 6,
             initial_phase=0):
        """
        """
        self.init_medium = mc[0]
        self.trident_angle = trident_angle
        self.theta_max = theta_max
        AA = 1 - np.cos(theta_max)
        self.nrays = nrays
        self.initial_phase = initial_phase
        vr = self.axial_ray.perpendicularDirection()
        vr = Vec3.rotate(vr, self.axial_ray.d, np.random.random() * np.pi * 2)

        bundle_dict = dict()
        total_tr_num = 0
        for i in range(self.nrays):
            # initialize n_rays number of random directed trident rays
            # theta = np.random.random() * self.theta_max
            theta = np.arccos(1 - AA * np.random.random())
            p1 = self.axial_ray.to_coordinate(np.cos(theta))

            # random angle on the ring by counter-clockwise rotation
            beta = np.random.random() * np.pi * 2
            v_ = Vec3.rotate(vr, self.axial_ray.d, beta)

            p_end = p1 + Vec3.normalize(v_) * np.sin(theta)
            pow_dire = p_end - self.center
            ray_helper = Ray(self.center, pow_dire)

            # two perpendicular directions for auxrays
            n_helper = ray_helper.perpendicularDirection()
            n2_helper = np.cross(pow_dire, n_helper)
            l = np.tan(self.trident_angle) * np.linalg.norm(pow_dire)
            assert l > 0
            a1_end = p_end + Vec3.normalize(n_helper) * l
            a2_end = p_end + Vec3.normalize(n2_helper) * l

            # calculate initial power at theta
            z = self.distance_z[LONGITUDINAL]
            pressure0 = self.ffa(z, theta)
            I0 = np.abs(pressure0)**2 / (2 * self.init_medium.Z[LONGITUDINAL])
            tr = Trident(
                self.center,
                pow_dire,
                self.center,
                a1_end - self.center,
                self.center,
                a2_end - self.center,
                I0,
                self.frequency,
                self.initial_phase,
                len0=z,
                el_id=self.el_id,
                ray_id=self.el_id * self.nrays + i,
                medium=mc[0],
                legacy=[],
                wave_type=LONGITUDINAL)
            power_limit = tr.P0 / 100  # TODO
            # https://docs.python.org/3/tutorial/datastructures.html
            tr_queue = deque()
            tr_queue.append(tr)
            total_tr_num += 1
            while len(tr_queue):
                # print("before reflect:", len(tr_queue))
                tnow = tr_queue.popleft()
                if tnow.P0 < power_limit:
                    # print("discarded")
                    continue
                if tnow.bundle_identifier not in bundle_dict:
                    bundle_dict[tnow.bundle_identifier] = []
                bundle_dict[tnow.bundle_identifier].append(tnow)
                
                tr_list = tnow.reflect(mc)
                for tr_ in tr_list:
                    tr_queue.append(tr_)
                    total_tr_num += 1
                tr_list2 = tnow.refract(mc)
                for tr_ in tr_list2:
                    tr_queue.append(tr_)
                    total_tr_num += 1
        print(f"total trident num: {total_tr_num}")
        return bundle_dict

    def Dfunc(self, theta):
        """
        helper function D used in `ffa` by convention
        f[x_] := BesselJ[1, k a Sin[x]]
        f'[x_] := 1/2 a k (BesselJ[0, a k Sin[x_]] - BesselJ[2, a k Sin[x_]]) Cos[x_]
        """
        if theta:
            var1 = self.ka[LONGITUDINAL] * np.sin(theta)
            return 2 * special.jv(1, var1) / var1
        else:
            return special.jv(0, 0) - special.jv(2, 0)

    def ffa(self, r, theta=None):
        """
        calculate pressure at point `r` using "Far Field Approximation"

        `p(r, theta) = area S0 exp(i k r)/(2 pi r) D(theta)`

        with: `area = pi Radius^2`
        """
        if theta is None:
            theta = Vec3.angle_between(r, self.axial_ray.d)
            r = np.linalg.norm(r)

        term2 = np.exp(1j * self.k[LONGITUDINAL] * r) / (2 * np.pi * r)
        c_pressure = self.area * self.S0 * term2 * self.Dfunc(theta)
        return c_pressure

    @cached_property
    def area(self):
        # area of flat transducer element
        return np.pi * self.radius**2

    @cached_property
    def k(self):
        # wave number, only used in transducer initialization
        return 2 * np.pi * self.frequency / self.init_medium.c

    @cached_property
    def ka(self):
        # ka = k * a
        return self.k * self.radius

    @cached_property
    def max_flux(self):
        # from Boris' code
        max_flux = integrate.quad(self.fluxfunc, 0, self.theta_max)
        return max_flux

    @cached_property
    def S0(self):
        """
        From Huub's ThreeRaysV2.m script
        computes S0 such that flat transducer element generates RequiredPower
        p(r, theta) = A S0 exp(i k r)/(2 pi r)  D(theta) with:
        A= area=pi Radius^2
        D(theta)= 2 J1(Radius*k sin(theta)) / (Radius*k*sin(theta))
        D = directivity = determines pressure field in direction theta
        J1=Besselfunction, n=1
        """
        Int = integrate.quad(self.fluxfunc, 0.00000001, np.pi / 2)
        S0 = np.sqrt(4 * np.pi * self.init_medium.Z[LONGITUDINAL] *
                     self.required_power / Int[0]) / self.area
        return S0

    @cached_property
    def wave_length(self):
        return self.init_medium.c / self.frequency

    @cached_property
    def distance_z(self):
        """
        < Biomedical Ultrasound > p158, needed to calculate initial power/intensity
        """
        return 2.3 * self.radius**2 / self.wave_length


class Transducer(list):
    """ HIFU Transducer: as list of elements """

    def __init__(self,
                 nature_focus=0,
                 actual_focus=0,
                 focus_diameter=0,
                 frequency=0,
                 element_coordinates=None,
                 element_radius=None,
                 element_power=None,
                 element_properties=None,
                 **kw):
        super().__init__()
        # self.sphere = Sphere(**sphere_configs)
        self.element_coordinates = element_coordinates
        self.nature_focus = nature_focus  # nature focus of the transducer sphere
        self.actual_focus = actual_focus  # ultrasound focus
        self.focus_diameter = focus_diameter  # diameter of the focus area
        self.element_radius = (element_radius or element_properties['radius']
                               )  # radius of transducer element
        self.element_power = element_power or element_properties['power']
        self.element_properties = element_properties
        self.frequency = frequency  # frequency of the US emitted

    def initialize(self,
                   init_medium,
                   n_rays=None,
                   trident_angle=None,
                   theta_max=None,
                   n_core=None,
                   verbose=False):
        """initialize transducer directly instead.

        initialize the transducer element and cast rays in init medium
        
        Arguments:
            init_medium {InitMedium} -- initial medium e.g. lossless and markoil
        
        Keyword Arguments:
            n_rays {int} -- number of ray per transducer (default: {None})
            trident_angle {float64} -- angle between pow_ray and aux_ray (default: {None})
            theta_max {float} -- angle between which rays are casted from transducer (default: {None})
            n_core {int} -- number of cpu core to use (default: {None})
            verbose {bool} -- whether to print info or not (default: {False})
        """

        self.init_medium = init_medium
        n_rays = n_rays
        trident_angle = trident_angle
        theta_max = theta_max
        initial_phase = 0
        seed = int(time.time())
        if verbose: print("random seed:", seed)
        # np.random.seed(18973894)
        np.random.seed(34839194)
        if n_core is not None:
            pool = Pool(n_core)
            async_results = []
            for i, co in enumerate(self.element_coordinates):
                ar = pool.apply_async(
                    te_init_wrapper,
                    args=(i, co, self.element_radius, self.element_power,
                          self.frequency, self.nature_focus, self.init_medium,
                          initial_phase, n_rays, trident_angle, theta_max,
                          verbose, True))
                async_results.append(ar)
            pool.close()
            pool.join()
            for ar in async_results:
                te = ar.get()
                self.append(te)
        else:
            for i, co in enumerate(self.element_coordinates):
                te = te_init_wrapper(
                    i, co, self.element_radius, self.element_power,
                    self.frequency, self.nature_focus, self.init_medium,
                    initial_phase, n_rays, trident_angle, theta_max, verbose)
                if verbose:
                    print("Initialized transducer #{}".format(te.el_id))
                self.append(te)

    def cast(self, mc=[]):
        """use cast in TElement class instead
        
        cast all the inital tridents towards a MediaComplex instance `mc`
        return dictionary of tridents sorted by bundle identifier string
        """

        bundle_dict = dict()
        for te in self:
            b = te.cast(mc)
            bundle_dict.update(b)
        return bundle_dict


def te_init_wrapper(el_id,
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
                    verbose,
                    parrallel=False):
    """multiprocessing should only be involved in HIFU module, use
    `pyHIFU.HIFU` instead
    
    Arguments:
        el_id {int} -- Element ID
        center {ndarray(1,3)} -- center of the transducer element
        radius {float} -- radius of transducer
        power {float} -- total power of transducer
        freq {float} -- frequence
        nature_f {ndarray(1,3)} -- nature focus
        init_medium {InitMedium} -- init medium
        initial_phase {float} -- init phase
        n_rays {int} -- number of rays from initial cast in init medium
        trident_angle {float} -- angle between pow_ray and aux_ray
        theta_max {float} -- angle between which rays are casted from transducer
        verbose {Bool} -- print info
    
    Keyword Arguments:
        parrallel {bool} -- running in multiprocessing (default: {False})
    
    Returns:
        TElement -- Transducer element instance
    """

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
    if verbose:
        pname = current_process().name
        print(f"TE #{te.el_id} initialized by {pname}")
    return te
