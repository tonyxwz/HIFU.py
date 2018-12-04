import abc
import numpy as np
from geometric.basic import Ray as GeoRay

        
class Ray(GeoRay):
    def __init__(self, start=None, direction=None):
        super().__init__(start, direction)

    @property
    def end(self):
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


class AuxRay(Ray):
    def __init__(self, start=None, direction=None,
            wave_type=None, polarization=None, initial_phase=None,
            medium=None):
        super().__init__(start, direction)
        self.wave_type = wave_type
        self.polarization = polarization
        self.initial_phase = initial_phase
        # self.end = 0 # calculate the end of the vector acc. medium
