import abc
import numpy as np
from sympy.geometry.line import LinearEntity3D, Ray3D
from sympy.geometry.line import Ray as SympyRay


class Ray(LinearEntity3D, SympyRay):
    """
    abstract class `segment`
    to be inherited by `PowRaySegment` and `AuxRaySegment`
    """
    def __init__(self, start=None, direction=None):
        direction = np.array(direction) 
        self.unit_vector = direction / np.linalg.norm(direction)
    
    def __new__(cls, start=None, direction=None, **kwargs):
        p1 = np.array(start)
        p2 = p1 + np.array(direction)
        return super().__new__(cls, p1, p2, **kwargs)

    def print_direction(self):
        print(self.direction)


class PowRay(Ray):
    """
    segments of a ray with power, a.k.a. main ray
    """
    def __init__(self, start=None, direction=None,
            wave_type=None, polarization=None, initial_phase=None,
            medium=None):
        super().__init__(start=start, direction=direction)
        self.wave_type = wave_type
        self.polarization = polarization
        self.initial_phase = initial_phase
        # self.end = 0 # calculate the end of the vector acc. medium

    @property
    def end(self):
        return 0


class AuxRay(Ray):
    """
    segment of an auxilliary ray
    """
    def __init__(self, start=None, direction=None,
            wave_type=None, polarization=None, initial_phase=None,
            medium=None):
        super().__init__(start=start, direction=direction)
        pass

    @property
    def end(self):
        return 0

