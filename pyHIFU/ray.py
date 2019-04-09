import abc
from copy import copy as shlwcopy
import numpy as np
from .geometric.lines import Ray as GeoRay
from .geometric.surfaces import Plane
from .physics.acoustics import *
# from pyHIFU.geometric.voxel import Voxel
from .physics import LONGITUDINAL, SHEAR, SOLID, LIQUID
from cached_property import cached_property
import warnings
warnings.simplefilter("error")

class Ray(GeoRay):
    def __init__(self,
                 start=None,
                 direction=None,
                 wave_type=LONGITUDINAL,
                 medium=None):
        super().__init__(start, direction)
        try:
            x_inv = 1 / self.d[0]
        except ZeroDivisionError:
            x_inv = np.inf

        try:
            y_inv = 1 / self.d[1]
        except ZeroDivisionError:
            y_inv = np.inf

        try:
            z_inv = 1 / self.d[2]
        except ZeroDivisionError:
            z_inv = np.inf

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
            if self.has_point(end):
                self._end = np.array(end)
                self._endt = (self._end[0] - self.start[0]) / self.d[0]
                self.terminated = True
            else:
                raise Exception("End point is not on ray")
        else:
            raise Exception("Ray is already terminated")

    @property
    def endt(self):
        """ end of ray in terms of distance travelled on the ray """
        if self.terminated:
            return self._endt
        else:
            raise Exception("Ray is not terminated")

    @endt.setter
    def endt(self, endt):
        if not self.terminated:
            self._end = self.to_coordinate(endt)
            self._endt = endt
            self.terminated = True
        else:
            raise Exception("Ray is already terminated")

    # @property
    # def wave_number(self):
    #     return
    @cached_property
    def constant1(self):
        pass

    @cached_property
    def constant2(self):
        pass

    @cached_property
    def constant3(self):
        pass

    @cached_property
    def constant4(self):
        pass

    @cached_property
    def constant5(self):
        pass

    @cached_property
    def constant6(self):
        pass

    @cached_property
    def constant7(self):
        pass

    # def FSolvePars(self):

    def v_reflect(self, interface):
        pass

    def v_refract(self, interface):
        pass


class PowRay(Ray):
    def __init__(self,
                 start=None,
                 direction=None,
                 frequency=0,
                 wave_type=LONGITUDINAL,
                 polarization=None,
                 initial_phase=None,
                 medium=None):
        super().__init__(
            start=start,
            direction=direction,
            wave_type=wave_type,
            medium=medium)
        self.frequency = frequency
        self.polarization = polarization
        self.initial_phase = initial_phase
        self.angular_frequency = 2 * np.pi * self.frequency
        self.wave_number = self.angular_frequency / self.medium.c[
            self.wave_type]
        self.wave_length = self.medium.c[self.wave_type] / self.frequency

    def snell(self):
        pass

    def get_phase_at(self, l):
        return (l / self.wave_length) * np.pi * 2

    @property
    def final_phase(self):
        return self.get_phase_at(self.endt)


class AuxRay(Ray):
    def __init__(self,
                 start=None,
                 direction=None,
                 wave_type=LONGITUDINAL,
                 medium=None):
        super().__init__(start, direction, wave_type=wave_type, medium=medium)
        # self.end = 0 # calculate the end of the vector acc. medium

    def snell(self):
        pass


class Trident(object):
    """
    power ray * 1 + auxiliary ray * 2
    """

    def __init__(self,
                 start_pow,
                 dire_pow,
                 start_aux1,
                 dire_aux1,
                 start_aux2,
                 dire_aux2,
                 I0,
                 frequency,
                 initial_phase,
                 len0=0,
                 el_id=None,
                 ray_id=None,
                 medium=None,
                 legacy=[],
                 wave_type=LONGITUDINAL):

        self.wave_type = wave_type
        self.frequency = frequency
        self.pow_ray = PowRay(
            start=start_pow,
            direction=dire_pow,
            frequency=frequency,
            initial_phase=initial_phase,
            wave_type=self.wave_type,
            medium=medium)
        self.aux_ray1 = AuxRay(
            start=start_aux1,
            direction=dire_aux1,
            wave_type=self.wave_type,
            medium=medium)
        self.aux_ray2 = AuxRay(
            start=start_aux2,
            direction=dire_aux2,
            wave_type=self.wave_type,
            medium=medium)
        self.medium = medium
        # self.med_id = med_id # medium index
        self.id = ray_id  # ray id
        self.el_id = el_id  # element index

        # Must copy to assign value, shallow copy is enough as str are inmutable
        self.history = shlwcopy(legacy)
        self.history.append(str(self.medium.id))

        self.I0 = I0  # initial intensity of pow_ray
        self.len0 = len0  # position (on pow_ray) at which initial intensity is calculate, aka distance_z
        self.P0 = self.I0 * self.get_area_at(len0)
        self.intersect_medium()
        # TODO
        self.PF = self.P0 * self.attenufactor(self.pow_ray.endt - self.len0)
        self.IF = self.get_intensity_at(self.pow_ray.endt)

    def intersect_medium(self):
        t = np.inf
        i = 0
        for i_, face in enumerate(self.medium.shape):
            t_ = face.intersect_line(line=self.pow_ray, require_t=True)
            if t_ is not None:
                if t_ < t and Vec3.greaterthan(t_, 0):
                    t = t_
                    i = i_
        self.pow_ray.endt = t
        self.exiting_face_index = i
        # TODO is face is a barrelshell, use the tangitial plane instead
        t1 = self.medium.shape[self.exiting_face_index].intersect_line(
            line=self.aux_ray1, require_t=True, as_plane=True)
        self.aux_ray1.endt = t1
        t2 = self.medium.shape[self.exiting_face_index].intersect_line(
            line=self.aux_ray2, require_t=True, as_plane=True)
        self.aux_ray2.endt = t2

    def next_medium(self, mc):
        """return the next medium

        Parameters
        ----------
        mc : MediaComplex
        
        Returns
        -------
        Medium
            next medium
        """
        for index, tf in enumerate(mc.adj_mtx[self.medium.id] == self.exiting_face_index):
            if tf:
                return mc[index]
        return None

    def reflect(self, mc):
        k = list()
        next_medium = self.next_medium(mc)
        # if there's no next medium, stop propagation
        if self.medium.is_init:
            # only consider refraction in init media
            return k
        if next_medium is None:
            # return k
            # TODO find better implementation
            next_medium = mc[self.medium.id - 1]
        T, R, Ph_refl, Ph_tr = coefficient_ll(self.pow_ray.d,
                                              self.medium.shape[self.exiting_face_index].n,
                                              self.medium.c[self.wave_type],
                                              next_medium.c[self.wave_type],
                                              self.medium.density,
                                              next_medium.density)
        assert Vec3.num_are_equal(T+R, 1)
        # reflected trident, same type
        pow_dir_1 = vray_reflected(self.pow_ray,
                                   self.medium.shape[self.exiting_face_index])
        aux1_dir_1 = vray_reflected(self.aux_ray1,
                                    self.medium.shape[self.exiting_face_index])
        aux2_dir_1 = vray_reflected(self.aux_ray2,
                                    self.medium.shape[self.exiting_face_index])

        # https://en.wikipedia.org/wiki/Reflection_phase_change
        if pow_dir_1 is not None and aux1_dir_1 is not None and aux2_dir_1 is not None:
            tr_1 = Trident(
                self.pow_ray.end,
                pow_dir_1,
                self.aux_ray1.end,
                aux1_dir_1,
                self.aux_ray2.end,
                aux2_dir_1,
                self.get_intensity_at(self.pow_ray.endt) * R,  # TODO coefficient
                self.frequency,
                self.get_phase_at(self.pow_ray.endt) + Ph_refl,  # TODO reflection shift
                el_id=self.el_id,
                ray_id=self.id,
                medium=self.medium,
                legacy=self.history,
                wave_type=self.wave_type)
            k.append(tr_1)

        # reflected trident, different type
        if self.medium.state == SOLID:
            new_wave_type = (self.pow_ray.wave_type + 1) % 2  # fancy
            pow_dir_2 = vray_reflected(
                self.pow_ray,
                self.medium[self.exiting_face_index],
                c1=self.medium.c[self.pow_ray.wave_type],
                c2=self.medium.c[new_wave_type])
            aux1_dir_2 = vray_reflected(
                self.aux_ray1,
                self.medium[self.exiting_face_index],
                c1=self.medium.c[self.aux_ray1.wave_type],
                c2=self.medium.c[new_wave_type])
            aux2_dir_2 = vray_reflected(
                self.aux_ray2,
                self.medium[self.exiting_face_index],
                c1=self.medium.c[self.aux_ray2.wave_type],
                c2=self.medium.c[new_wave_type])
            if pow_dir_2 is not None and aux1_dir_2 is not None and aux2_dir_2 is not None:
                tr_2 = Trident(
                    self.pow_ray.end,
                    pow_dir_2,
                    self.aux_ray1.end,
                    aux1_dir_2,
                    self.aux_ray2.end,
                    aux2_dir_2,
                    self.get_intensity_at(self.pow_ray.endt),  # TODO coefficient
                    self.frequency,
                    self.get_phase_at(self.pow_ray.endt),  # TODO reflection shift
                    el_id=self.el_id,
                    ray_id=self.id,
                    medium=self.medium,
                    legacy=self.history,
                    wave_type=new_wave_type)
                k.insert(tr_2.wave_type, tr_2)
        # tr_1 is longit. shear wave not possible in liquid, do nothing
        elif self.medium.state == LIQUID:
            tr_2 = None  # TODO remove not used

        return k

    def refract(self, mc):
        k = list()
        next_medium = self.next_medium(mc)
        # if there's no next medium, stop propagation
        if next_medium is None or next_medium.is_init:
            return k

        # reflected trident, same type
        pow_dir_1 = vray_refracted(self.pow_ray,
                                   self.medium.shape[self.exiting_face_index],
                                   self.medium.c[self.wave_type],
                                   next_medium.c[self.wave_type])
        aux1_dir_1 = vray_refracted(self.aux_ray1,
                                    self.medium.shape[self.exiting_face_index],
                                    self.medium.c[self.wave_type],
                                    next_medium.c[self.wave_type])
        aux2_dir_1 = vray_refracted(self.aux_ray2,
                                    self.medium.shape[self.exiting_face_index],
                                    self.medium.c[self.wave_type],
                                    next_medium.c[self.wave_type])

        # https://en.wikipedia.org/wiki/Reflection_phase_change
        if pow_dir_1 is not None and aux1_dir_1 is not None and aux2_dir_1 is not None:
            T, R, Ph_refl, Ph_tr = coefficient_ll(self.pow_ray.d,
                                                  self.medium.shape[self.exiting_face_index].n,
                                                  self.medium.c[self.wave_type],
                                                  next_medium.c[self.wave_type],
                                                  self.medium.density,
                                                  next_medium.density)
            assert Vec3.num_are_equal(T + R, 1)
            # total internal refraction
            tr_1 = Trident(
                self.pow_ray.end,
                pow_dir_1,
                self.aux_ray1.end,
                aux1_dir_1,
                self.aux_ray2.end,
                aux2_dir_1,
                self.get_intensity_at(self.pow_ray.endt) * T,
                self.frequency,
                self.get_phase_at(self.pow_ray.endt) + Ph_tr,
                el_id=self.el_id,
                ray_id=self.id,
                medium=next_medium,
                legacy=self.history,
                wave_type=self.wave_type)
            k.append(tr_1)

        # refracted trident, different type
        if next_medium.state == SOLID:
            new_wave_type = (self.pow_ray.wave_type + 1) % 2  # fancy
            pow_dir_2 = vray_refracted(
                self.pow_ray,
                self.medium[self.exiting_face_index],
                c1=self.medium.c[self.wave_type],
                c2=next_medium.c[new_wave_type])
            aux1_dir_2 = vray_refracted(
                self.aux_ray1,
                self.medium[self.exiting_face_index],
                c1=self.medium.c[self.wave_type],
                c2=next_medium.c[new_wave_type])
            aux2_dir_2 = vray_refracted(
                self.aux_ray2,
                self.medium[self.exiting_face_index],
                c1=self.medium.c[self.wave_type],
                c2=next_medium.c[new_wave_type])
            if pow_dir_2 is not None and aux1_dir_2 is not None and aux2_dir_2 is not None:
                # total internal refraction
                tr_2 = Trident(
                    self.pow_ray.end,
                    pow_dir_2,
                    self.aux_ray1.end,
                    aux1_dir_2,
                    self.aux_ray2.end,
                    aux2_dir_2,
                    self.get_intensity_at(self.pow_ray.endt),  # TODO coefficient
                    self.frequency,
                    self.get_phase_at(self.pow_ray.endt),  # TODO reflection shift
                    el_id=self.el_id,
                    ray_id=self.id,
                    medium=next_medium,
                    legacy=self.history,
                    wave_type=new_wave_type)
                k.insert(tr_2.wave_type, tr_2)
        # tr_1 is longit. shear wave not possible in liquid, do nothing
        elif self.medium.state == LIQUID:
            tr_2 = None
        return k

    @property
    def bundle_identifier(self):
        # TODO return unique bundle identifier according to properties
        transducer_string = 'tr_' + str(self.el_id)
        history_string = 'mh_' + '_'.join(self.history)
        type_string = "LONGIT" if self.wave_type == LONGITUDINAL else "SHEAR"

        return '_'.join([transducer_string, history_string, type_string])

    def __str__(self):
        return str(self.__dict__)

    def get_area_at_alt(self, distance):
        p0 = distance * self.pow_ray.unit_vector + self.pow_ray.p
        p1 = distance * self.aux_ray1.unit_vector + self.aux_ray1.p
        p2 = distance * self.aux_ray2.unit_vector + self.aux_ray2.p
        area = np.linalg.norm(np.cross(p1 - p0, p2 - p0)) / 2
        return area

    def get_area_at_heron(self, distance):
        """
        https://en.wikipedia.org/wiki/Heron%27s_formula
        """
        p0 = distance * self.pow_ray.unit_vector + self.pow_ray.p
        # new Plane instance too slow

        plane = Plane(p0, self.pow_ray.d)
        p1 = plane.intersect_line(self.aux_ray1)
        p2 = plane.intersect_line(self.aux_ray2)

        a = np.linalg.norm(p0 - p1)
        b = np.linalg.norm(p1 - p2)
        c = np.linalg.norm(p2 - p0)

        s = (a + b + c) / 2
        area = np.sqrt(s*(s-a)*(s-b)*(s-c))
        return area
    
    def get_area_at(self, distance):
        """
        https://proofwiki.org/wiki/Norm_of_Vector_Cross_Product
        """
        p0 = distance * self.pow_ray.unit_vector + self.pow_ray.p
        # new Plane instance too slow

        plane = Plane(p0, self.pow_ray.d)
        p1 = plane.intersect_line(self.aux_ray1)
        p2 = plane.intersect_line(self.aux_ray2)

        try:
            area = np.linalg.norm(np.cross(p1 - p0, p2 - p0)) / 2
        except RuntimeWarning:
            print(p0, p1, p2)
        return area

    def get_intensity_at(self, d):
        """
        `d`: distance (scalar)
        """
        A1 = self.get_area_at(d)
        # TODO add attenuation here (or rename to `get_intensity_at`)
        I1 = self.P0 * self.attenufactor(d - self.len0) / A1
        # Q = I1 * 2 * attenuation
        return I1

    def get_phase_at(self, d):
        return self.pow_ray.get_phase_at(d)

    def attenufactor(self, s):
        return np.exp(-2 * self.medium.attenuation[self.wave_type] * s)
