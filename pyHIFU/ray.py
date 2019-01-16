import abc

import numpy as np

from .geometric.lines import Ray as GeoRay
from .geometric.surfaces import Plane
from pyHIFU.physics import snell, fresnel
from pyHIFU.geometric.voxel import Voxel

class Ray(GeoRay):
    def __init__(self, start=None, direction=None):
        super().__init__(start, direction)

    @property
    def end(self):
        pass

    def v_reflect(self, interface):
        pass

    def v_refract(self, interface):
        pass


class PowRay(Ray):
    def __init__(self, start=None, direction=None,
                 wave_type=None, polarization=None, initial_phase=None,
                 medium=None):
        super().__init__(start, direction)
        self.wave_type = wave_type
        self.polarization = polarization
        self.initial_phase = initial_phase
        # self.end = 0 # calculate the end of the vector acc. medium

    def powsnell(self):
        pass


class AuxRay(Ray):
    def __init__(self, start=None, direction=None,
                 wave_type=None, polarization=None, initial_phase=None,
                 medium=None):
        super().__init__(start, direction)
        self.wave_type = wave_type
        self.polarization = polarization
        self.initial_phase = initial_phase
        # self.end = 0 # calculate the end of the vector acc. medium

    def auxsnell(self):
        pass

class Trident(object):
    """
    power ray * 1 + auxilliary ray * 2
    """

    def __init__(self, start_pow, dire_pow,
                 start_aux1, dire_aux1,
                 start_aux2, dire_aux2,
                 el_idx=None, ray_idx=None, med_idx=None):
            # medium)
        self.pow_ray = PowRay(start=start_pow, direction=dire_pow)
        self.aux_ray1 = AuxRay(start=start_aux1, direction=dire_aux1)
        self.aux_ray2 = AuxRay(start=start_aux2, direction=dire_aux2)
        self.el_idx=el_idx # element index
        self.id = ray_idx # ray id
        self.med_idx = med_idx # medium index

    def __str__(self):
        return str(self.__dict__)

    def get_area_at(self, distance):
        p0 = distance * self.pow_ray.unit_vector + self.pow_ray.p
        # new Plane instance too slow

        plane = Plane(p0, self.pow_ray.d)
        p1 = plane.intersect_line(self.aux_ray1)
        p2 = plane.intersect_line(self.aux_ray2)

        return np.linalg.norm(np.cross(p1-p0, p2-p0)) / 2

    def get_power_at(self, distance):
        pass

    def trisnell(self):
        # get "snelled" tridents, call `powsnell` and `auxsnell`
        pass

    def intersect_voxel(self, v:Voxel):
        pass