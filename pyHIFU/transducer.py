from collections import deque
import time
import numpy as np
from .geometric.surfaces import Sphere
from .geometric.lines import Ray
from .geometric.surfaces import Plane
from .geometric.curves import Ring
from .ray import Trident
from .io.config import readjson
from .geometric.vec3 import Vec3
from scipy import special, integrate
from pyHIFU.physics import LONGITUDINAL


class TElement(list):
    """ tranducer element class """
    def __init__(self, el_id, center, radius=0, power=0, freq=0,
                 nature_f=np.array([0,0,0]), **kw):
        """ only assgin necessary parameters here """
        super().__init__()  # len(self) == 0
        self.el_id = el_id
        self.center = center  # also the start of all tridents
        self.radius = radius
        self.required_power = power  # needed in ffa
        self.frequency = freq
        self.nature_f = nature_f
        self.axial_ray = Ray(self.center, self.nature_f-self.center)

    def initialize(self, init_medium, initial_phase, n=100, trident_angle=1e-4, theta_max=np.pi/6):
        """
        initialize all trident rays until they hit markoil interface
        `init_medium`: e.g. lossless / markoil
        `n`: number of rays per transducer
        `trident_angle`: the angle between powray and auxray
        `sampling_angle`: theta_max
        """
        self.trident_angle = trident_angle
        self.theta_max = theta_max
        self.n_rays = n
        self.init_medium = init_medium
        self.fluxfunc = lambda x: np.sin(x) * (2 * special.jv(1, self.ka * x) / (self.ka * x) )**2
        self.initial_phase = initial_phase*np.random.random()
        vr = self.axial_ray.perpendicularDirection()
        vr = Vec3.rotate(vr, self.axial_ray.d, np.random.random()*np.pi*2)
        for i in range(self.n_rays):
            # initialize n_rays number of random directed trident rays
            theta = np.random.random() * self.theta_max
            p1 = self.axial_ray.to_coordinate(np.cos(theta))

            # random angle on the ring by counter-clockwise rotation
            beta = np.random.random() * np.pi * 2
            v_ = Vec3.rotate(vr, self.axial_ray.d, beta)

            p_end = p1 + Vec3.normalize(v_)*np.sin(theta)
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
            I0 = np.abs(pressure0)**2 / (2*self.init_medium.Z[LONGITUDINAL])
            self.append(Trident(self.center, pow_dire,
                                self.center, a1_end-self.center,
                                self.center, a2_end-self.center,
                                I0, self.frequency, self.initial_phase, len0=z,
                                el_id=self.el_id, ray_id=self.el_id*self.n_rays+i,
                                medium=init_medium, legacy=[],
                                wave_type=LONGITUDINAL))

    # .-----------------------.
    # | short-hand properties |
    # '-----------------------'
    @property
    def area(self):
        # area of flat transducer element
        return np.pi * self.radius**2
    @property
    def k(self):
        # wave number, only used in transducer initialization
        return 2*np.pi*self.frequency / self.init_medium.c
    @property
    def ka(self):
        # ka = k * a
        return self.k * self.radius
    @property
    def max_flux(self):
        # from Boris' code
        max_flux = integrate.quad(self.fluxfunc, 0, self.theta_max)
        return max_flux
    @property
    def S0(self):
        """
        From Huub's ThreeRaysV2.m script
        computes S0 such that flat transducer element generates RequiredPower
        p(r, theta) = A S0 exp(i k r)/(2 pi r)  D(theta) with:
        A= area=pi Radius^2
        D(theta)= 2 J1(Radius*k sin(theta)) / (Radius*k*sin(theta))
        D = directivity = determines pressure field in direction theta
        J1=Besselfunction, n=1
        TODO find which is CORRECT formula
        """
        Int = integrate.quad(self.fluxfunc, 0.00000001, np.pi/2)
        S0 = np.sqrt(4*np.pi*self.init_medium.Z[LONGITUDINAL]*self.required_power/Int[0])/self.area
        return S0
    @property
    def wave_length(self):
        return self.init_medium.c / self.frequency
    @property
    def distance_z(self):
        """
        < Biomedical Ultrasound > p158, needed to calculate initial power/intensity
        """
        return 2.3 * self.radius**2 / self.wave_length
    # .----------------------.
    # |   end of properties  |
    # '----------------------'
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


class Transducer(list):
    """ HIFU Transducer: as list of elements """
    def __init__(self,
                 nature_focus=0, actual_focus=0, focus_diameter=0, frequency=0,
                 element_properties=None, element_coordinates=None,
                 element_init_paras=None):
        super().__init__()
        # self.sphere = Sphere(**sphere_configs)
        self.element_coordinates = element_coordinates
        self.nature_focus = nature_focus              # nature focus of the transducer sphere
        self.actual_focus = actual_focus              # ultrasound focus
        self.focus_diameter = focus_diameter          # diameter of the focus area
        self.element_properties = element_properties  # radius of transducer element
        self.frequency = frequency                    # frequency of the US emitted

        for i,co in enumerate(self.element_coordinates):
            self.append(TElement(i, co,
                                 radius=self.element_properties["radius"],
                                 power=self.element_properties["power"],
                                 freq=self.frequency,
                                 nature_f=self.nature_focus))

    def initialize(self, init_medium, n_rays=None, trident_angle=None, theta_max=None):
        self.init_medium = init_medium
        n_rays = n_rays
        trident_angle = trident_angle
        theta_max = theta_max

        interface = self.init_medium.shape[0]
        seed = int(time.time())
        print("random seed:", seed)
        np.random.random(seed)
        for te in self:
            te.initialize(self.init_medium, initial_phase=0,
                          n=n_rays,
                          trident_angle=trident_angle,
                          theta_max=theta_max)

            for tr in te:
                # TODO move set end point to casting
                tr.pow_ray.end = interface.intersect_line(tr.pow_ray)
                # tr.aux_ray1.end = interface.intersect_line(tr.aux_ray1)
                # tr.aux_ray2.end = interface.intersect_line(tr.aux_ray2)

    def cast(self, mc=[]):
        """ cast all the inital tridents towards a MediaComplex instance `mc`
        return dictionary of tridents sorted by bundle identifier string
        """

        bundle_dict = dict()
        for te in self:
            if len(te) == 0:
                raise Exception("Must initialize Transducer to cast rays")
            for tr in te:
                # https://docs.python.org/3/tutorial/datastructures.html#using-lists-as-stacks
                tr_stack = [tr]
                while len(tr_stack):
                    tnow = tr_stack.pop()
                    if not tnow.bundle_identifier in bundle_dict:
                        bundle_dict[tnow.bundle_identifier] = []
                    bundle_dict[tnow.bundle_identifier].append(tr)
                    # TODO
                    # t1, t2 = tnow.reflect(mc)
                    # t3, t4 = tnow.refract(mc)
                    # t1.end = ... t2.end = ... t3.end = ... t4.end = ...
                    # tr_stack.append(t1)
                    # tr_stack.append(t2)
                    # tr_stack.append(t3)
                    # tr_stack.append(t4)
        return bundle_dict
