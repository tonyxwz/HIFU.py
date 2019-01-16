import numpy as np
from .geometric.surfaces import Sphere
from .geometric.lines import Ray
from .geometric.surfaces import Plane
from .geometric.curves import Ring
from .ray import Trident
from .io.config import readjson
from .geometric.vec3 import Vec3


class TElement(list):
    """ tranducer element class """
    def __init__(self, center, el_idx, trident_angle=1e-4,
                 n_rays=500, sampling_angle=np.pi/4, **kw):
        super().__init__()
        self.center = center # also the start of all tridents
        self.el_idx = el_idx
        self.sampling_sphere = Sphere(center, 1, np.array([0,0,0])-self.center,
                                       angle=np.pi/2)
        self.axial_ray = Ray(center, np.array([0,0,0])-self.center)
        v_radial_origin = self.axial_ray.perpendicularDirection()
        for i in range(n_rays):
            alpha = np.random.random() * sampling_angle
            p1 = self.axial_ray.to_coordinate(np.cos(alpha))
            # plane1 = Plane(p1, axial_ray.d)
            # ring = Ring(p1, np.sin(alpha), axial_ray.d)
            # random angle on the ring by counter-clockwise rotation
            beta = np.random.random() * np.pi * 2
            v_ = Vec3.rotate(v_radial_origin, self.axial_ray.d, beta)

            p_end = p1 + Vec3.normalize(v_)*np.sin(alpha)
            pow_dire = p_end - self.center
            ray_helper = Ray(self.center, pow_dire)
            # two perpendicular directions for auxrays
            n_helper = ray_helper.perpendicularDirection()
            n2_helper = np.cross(pow_dire, n_helper)
            l = np.tan(trident_angle) * np.linalg.norm(pow_dire)
            assert l > 0
            a1_end = p_end + Vec3.normalize(n_helper) * l
            a2_end = p_end + Vec3.normalize(n2_helper) * l
            
            self.append(Trident(self.center, pow_dire,
                                self.center, a1_end-self.center,
                                self.center, a2_end-self.center,
                                el_idx=self.el_idx, ray_idx=self.el_idx*n_rays+i,
                                med_idx=0)) # med_idx = 0: lossless
            
    @staticmethod
    def initial_trident(start, angle):
        pass


    def cast(self, n):
        pass
        

class Transducer(list):
    def __init__(self, config):
        super().__init__()
        self.sphere = Sphere(**config['sphere'])
        el_list = config['element_list']
        for i,el in enumerate(el_list):
            # TODO: set ray number
            self.append(TElement(el["coordinates"], i))

    def initialize(self):
        pass

