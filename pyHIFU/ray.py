import abc
import numpy as np
from .geometric.lines import Ray as GeoRay
from .geometric.surfaces import Plane
from pyHIFU.physics.acoustics import snell#, fresnel
# from pyHIFU.geometric.voxel import Voxel
from pyHIFU.physics import LONGITUDINAL, SHEAR, SOLID, LIQUID

class Ray(GeoRay):
    def __init__(self, start=None, direction=None, wave_type=LONGITUDINAL, medium=None):
        super().__init__(start, direction)
        
        self.wave_type = wave_type
        
        self.medium = medium  # medium instance
        self.terminated = False
    
    def find_terminal(self, mc):
        """ 
        `mc`: MediaComplex
        """
        pass

    def set_end(self, end):
        if not self.terminated:
            self.end = np.array(end)
            self.terminated = True
        else:
            raise Exception("Ray is already terminated")


    # @property
    # def wave_number(self):
    #     return 
    @property
    def constant1(self):
        pass
    @property
    def constant2(self):
        pass
    @property
    def constant3(self):
        pass
    @property
    def constant4(self):
        pass
    @property
    def constant5(self):
        pass
    @property
    def constant6(self):
        pass
    @property
    def constant7(self):
        pass

    # def FSolvePars(self):
    
    def v_reflect(self, interface):
        pass

    def v_refract(self, interface):
        pass


class PowRay(Ray):
    def __init__(self, start=None, direction=None, frequency=0,
                 wave_type=LONGITUDINAL, polarization=None, initial_phase=None,
                 medium=None):
        super().__init__(start=start, direction=direction,
                         wave_type=wave_type, medium=medium)
        self.frequency = frequency
        self.polarization = polarization
        self.initial_phase = initial_phase
        self.angular_frequency = 2 * np.pi * self.frequency
        self.wave_number = self.angular_frequency / self.medium.c[self.wave_type]
        # self.end = 0 # calculate the end of the vector acc. medium

    def snell(self):
        pass


class AuxRay(Ray):
    def __init__(self, start=None, direction=None,
                 wave_type=LONGITUDINAL, medium=None):
        super().__init__(start, direction, wave_type=wave_type, medium=medium)
        # self.end = 0 # calculate the end of the vector acc. medium

    def snell(self):
        pass


class Trident(object):
    """
    power ray * 1 + auxilliary ray * 2
    """

    def __init__(self, start_pow, dire_pow,
                 start_aux1, dire_aux1,
                 start_aux2, dire_aux2,
                 I0, len0, frequency, initial_phase,
                 el_id=None, ray_id=None, medium=None,
                 wave_type=LONGITUDINAL):
            # medium)
        self.wave_type = wave_type

        self.pow_ray = PowRay(start=start_pow, direction=dire_pow,
                              frequency=frequency, initial_phase=initial_phase,
                              wave_type=self.wave_type, medium=medium)
        self.aux_ray1 = AuxRay(start=start_aux1, direction=dire_aux1,
                               wave_type=self.wave_type, medium=medium)
        self.aux_ray2 = AuxRay(start=start_aux2, direction=dire_aux2,
                               wave_type=self.wave_type, medium=medium)
        self.medium = medium
        self.id = ray_id # ray id
        self.el_id=el_id # element index
        # self.med_id = med_id # medium index
        self.history = list()

        self.I0 = I0  # initial intensity of pow_ray
        self.len0 = len0  # position (on pow_ray) at which initial intensity is calculate
        self.P0 = self.I0 * self.get_area_at(len0)

    def __str__(self):
        return str(self.__dict__)

    def get_area_at(self, distance):
        """
        https://proofwiki.org/wiki/Norm_of_Vector_Cross_Product
        """
        p0 = distance * self.pow_ray.unit_vector + self.pow_ray.p
        # new Plane instance too slow

        plane = Plane(p0, self.pow_ray.d)
        p1 = plane.intersect_line(self.aux_ray1)
        p2 = plane.intersect_line(self.aux_ray2)

        area = np.linalg.norm(np.cross(p1-p0, p2-p0)) / 2
        return area

    def get_intensity_at(self, d):
        """
        `d`: distance (scalar)
        """
        A1 = self.get_area_at(d)
        # TODO add attenuation here (or rename to `get_intensity_at`)
        I1 = self.P0 * self.attenufactor(d) / A1
        # Q = I1 * 2 * attenuation
        return I1

    def attenufactor(self, s):
        return np.exp(-2*self.medium.attenuation[self.wave_type]*s)

    def trisnell(self, boundary):
        # get "snelled" tridents, call `powsnell` and `auxsnell`
        pass
