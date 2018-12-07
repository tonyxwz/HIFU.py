import numpy as np
from pyHIFU.ray import PowRay, AuxRay
from pyHIFU.geometric.basic import Plane


class Trident(object):
    """
    power ray * 1 + auxilliary ray * 2
    """
    def __init__(self, start_pow, dire_pow,
            start_aux1, dire_aux1,
            start_aux2, dire_aux2,
            tranducer_el_no=None):
            # medium)
        self.pow_ray = PowRay(start=start_pow, direction=dire_pow)
        self.aux_ray1 = AuxRay(start=start_aux1, direction=dire_aux1)
        self.aux_ray2 = AuxRay(start=start_aux2, direction=dire_aux2)

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


# def vector_