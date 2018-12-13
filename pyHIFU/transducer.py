import numpy as np
from .geometric.surfaces import Sphere
from .ray import Trident
from .io.config import readjson


class TElement(list):
    """ tranducer element class """
    def __init__(self, coor, n_rays, el_idx, **kw):
        super().__init__()
        self.coor = coor
        for i in range(n_rays):
            self.append(Trident(coor, [0,1,0],
                                coor, [0,1,1e-4],
                                coor, [1e-4,1,0],
                                el_idx=el_idx, ray_idx=i))

    def cast(self, n):
        pass
        

class Transducer(list):
    def __init__(self, config):
        super().__init__()
        self.sphere = Sphere(**config['sphere'])
        el_list = config['element_list']
        for i,el in enumerate(el_list):
            # TODO: set ray number
            self.append(TElement(el["coordinates"], 1, i))

    def initialize(self):
        pass

