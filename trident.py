import numpy as np
from sympy import Plane, Line, Ray
from ray import PowRay, AuxRay
# from geometric.basic import Vector, Plane, Line


class Trident(object):
    """
    power ray * 1 + auxilliary ray * 2
    """
    def __init__(self, start_pow, dire_pow,
            start_aux1, dire_aux1,
            start_aux2, dire_aux2):
            # medium)
        self.pow_ray = PowRay(start=start_pow, direction=dire_pow)
        self.aux_ray1 = AuxRay(start=start_aux1, direction=dire_aux1)
        self.aux_ray2 = AuxRay(start=start_aux2, direction=dire_aux2)

    def get_area_at(self, distance):
        p0 = distance * self.pow_ray.unit_vector + np.array(self.pow_ray.source, dtype=np.float)
        # new Plane instance too slow
        plane = Plane(p0, normal_vector=self.pow_ray.direction_ratio)
        p1 = np.array(plane.intersection(self.aux_ray1)[0], dtype=np.float)
        p2 = np.array(plane.intersection(self.aux_ray2)[0], dtype=np.float)
        
        return np.linalg.norm(np.cross(p1-p0, p2-p0))
        
    def get_power_at(self, distance):
        pass


# def vector_