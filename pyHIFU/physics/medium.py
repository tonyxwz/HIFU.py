import numpy as np
from pyHIFU.io.config import readjson
from pyHIFU.geometric.volume import Ball, Cuboid, Cylinder

SHAPE_FUNC_DICT = {'ball': Ball, 'cuboid': Cuboid, 'cylinder': Cylinder}


class Medium(object):
    """ Medium """
    """ def __init__(self, type=None, name=None,
            speed=None, density=None, attenuation=None, absorption=None,
            frequency=None, faces=None, adj_mtx=None, config_file=None,
            **kwargs): """

    def __init__(self, phyinfo=None, geoinfo=None, **kw):
        # TODO: initialize physics info
        shape_type = geoinfo['shape_type']
        # print("Medium::geoinfo[]:", geoinfo['parameters'])
        self.shape = SHAPE_FUNC_DICT[shape_type](**geoinfo['parameters'])
        self.index = kw['index']


class MediaComplex(list):
    """
    a list of all the media that is present in the HIFU system
    use adjacency matrix to find neighbours for acoustics related calculation
    """

    # def __init__(self, media_list, adj_mtx):
    def __init__(self, config_file_path):
        super().__init__()
        config_json = readjson(json_path=config_file_path)
        self.__build(config_json["medium_list"])
        # self.adj_mtx = np.eye(len(self))

    def __build(self, med_list_json):
        for item in med_list_json:
            phyinfo = item["physics"]
            geoinfo = item["geometry"]
            n = len(self)
            self.append(Medium(phyinfo=phyinfo, geoinfo=geoinfo, index=n))
        
        # adj_mtx[i,j] = [1,2] : 1,2 are the indices of the faces of self[i]
        # adjacent to self[j]
        self.adj_mtx = list()
        
        # now you have the list ready, time to build adjacency matrix
        # traverse all media (double compare loop)
        for i,item1 in enumerate(self):
            l = list()
            for j,item2 in enumerate(self):
                if i == j:
                    l.append([])
                else:
                    l.append(item1.shape.adj_at(item2.shape))
            self.adj_mtx.append(l)

    def find_next(self, med_idx, side_idx):
        """
        find the index of next medium, given the index of current medium
        and outgoing side `side`
        """
        # TODO
        r = []
        for i,v in enumerate(self.adj_mtx[med_idx]):
            if side_idx in v:
                r.append(i)
        return r


if __name__ == "__main__":
    mc = MediaComplex('data/example2.json')
    print(mc.find_next(1, 4))