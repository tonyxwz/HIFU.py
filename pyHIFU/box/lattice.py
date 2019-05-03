import numpy as np


class Lattice:
    def __init__(self, p, l):
        self.o1 = np.array(p)
        self.a = np.array([l[0], 0, 0])
        self.b = np.array([0, l[1], 0])
        self.c = np.array([0, 0, l[2]])
        self.center = self.o1 + (self.a+self.b+self.c) / 2  # for calculating intensity
        self.o2 = self.o1 + self.a + self.b + self.c

    def rayBoxIntersection(self, ray):
        """
        TODO move this method to `Cuboid` class to avoid duplicated definition
        find start (tmin) and end (tmax) distances of the ray in the lattrix system
        `ray`: must be a instance of `pyHIFU.ray.Ray` (PowRay / AuxRay)
        return in and out points

        This program routine will be called on every ray in the system thus needs
        optimizing. (method only hold for AABBs "Axis Aligned Bounding Boxes")

        1. http://www.cs.utah.edu/~awilliam/box/box.pdf

        2. https://tavianator.com/fast-branchless-raybounding-box-intersections/
        """
        # vectorized operation
        txyz1 = (self.o1 - ray.p) * ray.d_inv
        txyz2 = (self.o2 - ray.p) * ray.d_inv

        tmins = np.min([txyz1, txyz2], axis=0)
        tmaxs = np.max([txyz1, txyz2], axis=0)

        tmin = np.max(tmins)  # max of the 3 mins
        tmax = np.min(tmaxs)  # min of the 3 maxs

        if tmin < tmax:
            if tmin > ray.endt or tmax < 0:
                return None, None
            else:
                if tmin < 0:
                    tmin = 0
                if tmax > ray.endt:
                    tmax = ray.endt
                return tmin, tmax
        else:
            return None, None

