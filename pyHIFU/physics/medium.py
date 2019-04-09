"""
Material only contains physics properties,
Medium inherit from Material and have shape info
"""
import numpy as np
from cached_property import cached_property

from pyHIFU.geometric.surfaces import Plane, Sphere
from pyHIFU.geometric.volumes import Ball, Cuboid, Cylinder
from pyHIFU.io.config import readjson
from pyHIFU.physics import LIQUID, LONGITUDINAL, SHEAR, SOLID

from pyHIFU.physics.material import Material

SHAPE_CLASS_DICT = {
    'ball': Ball,
    'cuboid': Cuboid,
    'cylinder': Cylinder,
    'plane': Plane
}

BOUNDARY_CLASS_DICT = {'plane': Plane, 'sphere': Sphere}

markoil_properties = {
    'material_name': 'markoil',
    'state': LIQUID,
    'density': 1070.0,
    'cL': 1430.0,
    'absorption': 1.04,
    'attenuationL': 1.04,
    'heat_capacity': 4200.0,
    'thermal_conductivity': 0.5
}

lossless_properties = {
    'material_name': 'lossless',
    'state': LIQUID,
    'cL': 1380.0,
    'density': 1030.0,
    'attenuationL': 0,
    'absorption': 0
}

air_properties = {

}


class Medium(Material):
    def __init__(self,
                 material_name=None,
                 med_name=None,
                 med_id=None,
                 geometry=None,
                 is_init=False,
                 **kw):

        super().__init__(material_name=material_name, **kw)
        self.id = med_id
        self.name = med_name
        self.is_init = is_init

        shape_type = geometry['shape_type']
        self.shape = SHAPE_CLASS_DICT[shape_type](**geometry['parameters'])


class MediaComplex(list):
    """
    a list of all the media that is present in the HIFU system
    use adjacency matrix to find neighbours for acoustics related calculation
    """

    def __init__(self, medium_list=None, config_file_path=None):
        super().__init__()
        if medium_list is None:
            medium_list = readjson(json_path=config_file_path)["medium_list"]

        for item in medium_list:
            if 'is_init' in item and item['is_init']:
                for m in self:
                    m.id += 1
                self.insert(0, Medium(**item, med_id=0))
            else:
                n = len(self)
                self.append(Medium(**item, med_id=n))

        # adj_mtx[i,j] = [1,2] : 1,2 are the indices of the faces of self[i]
        # adjacent to self[j]
        self.adj_mtx = np.ones((len(self), len(self))) * -1

        # traverse all media (double compare loop)
        for i, item1 in enumerate(self):
            for j, item2 in enumerate(self):
                if i != j:
                    self.adj_mtx[i, j] = item1.shape.adj_at(item2.shape)

    @staticmethod
    def from_config(config):
        pass

    def find_next(self, med_id, side_id):
        """
        find the index of next medium, given the index of current medium
        and outgoing side `side`
        """
        # TODO: if there's no adjacency, it is adjacent to the air.
        r = []
        for i, v in enumerate(self.adj_mtx[med_id]):
            if side_id in v:
                r.append(i)
        return r

    def find_next_init(self):
        """ not necessary any more """
        for i, v in enumerate(self):
            for j, f in enumerate(v):
                if self[0].boundary == f:
                    return i, j


if __name__ == "__main__":
    mc = MediaComplex(config_file_path='../../data/case2.json')
    print(mc.adj_mtx)
