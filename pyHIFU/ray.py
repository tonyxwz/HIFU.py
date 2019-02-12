import abc
from copy import copy as shlwcopy
import numpy as np
from .geometric.lines import Ray as GeoRay
from .geometric.surfaces import Plane
from pyHIFU.physics.acoustics import snell#, fresnel
# from pyHIFU.geometric.voxel import Voxel
from pyHIFU.physics import LONGITUDINAL, SHEAR, SOLID, LIQUID

class Ray(GeoRay):
    def __init__(self, start=None, direction=None, wave_type=LONGITUDINAL, medium=None):
        super().__init__(start, direction)
        try: x_inv = 1 / self.d[0]
        except ZeroDivisionError: x_inv = np.inf

        try: y_inv = 1 / self.d[1]
        except ZeroDivisionError: y_inv = np.inf

        try: z_inv = 1 / self.d[2]
        except ZeroDivisionError: z_inv = np.inf

        self.d_inv = np.array([x_inv, y_inv, z_inv])

        self.wave_type = wave_type

        self.medium = medium  # medium instance
        self.terminated = False

    def find_terminal(self, mc):
        """
        `mc`: MediaComplex
        """
        pass

    @property
    def end(self):
        if self.terminated:
            return self._end
        else:
            raise Exception("Ray is not terminated")
    @end.setter
    def end(self, end):
        if not self.terminated:
            self._end = np.array(end)
            self.terminated = True
        else:
            raise Exception("Ray is already terminated")
    @property
    def endt(self):
        """ end of ray in terms of distance travelled on the ray """
        if self.terminated:
            return (self.end[0]-self.start[0]) / self.d[0]
        else:
            raise Exception("Ray is not terminated")
    @endt.setter
    def endt(self, endt):
        if not self.terminated:
            self._end = self.to_coordinate(endt)
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
        self.wave_length = self.medium.c[self.wave_type] / self.frequency

    def snell(self):
        pass

    def get_phase_at(self, l):
        return (l / self.wave_length) * np.pi * 2

    @property
    def final_phase(self):
        return self.get_phase_at(self.endt)


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
                 I0, frequency, initial_phase, len0=0,
                 el_id=None, ray_id=None, medium=None, legacy=[],
                 wave_type=LONGITUDINAL):

        self.wave_type = wave_type

        self.pow_ray = PowRay(start=start_pow, direction=dire_pow,
                              frequency=frequency, initial_phase=initial_phase,
                              wave_type=self.wave_type, medium=medium)
        self.aux_ray1 = AuxRay(start=start_aux1, direction=dire_aux1,
                               wave_type=self.wave_type, medium=medium)
        self.aux_ray2 = AuxRay(start=start_aux2, direction=dire_aux2,
                               wave_type=self.wave_type, medium=medium)
        self.medium = medium
        # self.med_id = med_id # medium index
        self.id = ray_id # ray id
        self.el_id=el_id # element index

        # Must copy to assign value, shallow copy is enough as str are inmutable
        self.history = shlwcopy(legacy)
        self.history.append(str(self.medium.id))

        self.I0 = I0  # initial intensity of pow_ray
        self.len0 = len0  # position (on pow_ray) at which initial intensity is calculate, aka distance_z
        self.P0 = self.I0 * self.get_area_at(len0)

    @property
    def bundle_identifier(self):
        # TODO return unique bundle identifier according to properties
        transducer_string = 'tr_' + str(self.el_id)
        history_string = 'mh_' + '_'.join(self.history)
        type_string = "LONGIT" if self.wave_type == LONGITUDINAL else "SHEAR"

        return '_'.join([transducer_string, history_string, type_string])

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
        I1 = self.P0 * self.attenufactor(d-self.len0) / A1
        # Q = I1 * 2 * attenuation
        return I1

    def get_phase_at(self, d):
        return self.pow_ray.get_phase_at(d)

    def attenufactor(self, s):
        return np.exp(-2*self.medium.attenuation[self.wave_type]*s)

    def trisnell(self, boundary):
        # get "snelled" tridents, call `powsnell` and `auxsnell`
        pass
