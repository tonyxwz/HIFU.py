from pyHIFU.geometric.volumes import Cuboid
from pyHIFU.ray import Trident
import numpy as np


class Box(Cuboid):
    def __init__(self, x1, y1, z1, x2, y2, z2, l, n_trd=256):
        """
        `[xmin, ymin, zmin]` and `[xmax, ymax, zmax]` defines the diagonal of the Box

        `l`: the edge length of lattice
        """
        xmin = min(x1, x2); xmax = max(x1, x2)
        ymin = min(y1, y2); ymax = max(y1, y2)
        zmin = min(z1, z2); zmax = max(z1, z2)

        o = [xmin, ymin, zmin]
        a = [xmax-xmin, 0, 0]
        b = [0, ymax-ymin, 0]
        c = [0, 0, zmax-zmin]
        super().__init__(o, a, b, c)  # self.o1, self.o2

        self.l = l  # lattice step
        self.n_trd = n_trd  # num of transducers
        nx = np.int(np.ceil((xmax-xmin) / self.l))
        ny = np.int(np.ceil((ymax-ymin) / self.l))
        nz = np.int(np.ceil((zmax-zmin) / self.l))
        self.lattrix = list()  # LAT(tice ma)TRIX 3D
        for ix in range(nx):
            x = xmin + ix * self.l
            P = list()  # 2D (Plane)
            for iy in range(ny):
                y = ymin + iy * self.l
                R = list()  # 1D (Row)
                for iz in range(nz):
                    z = zmin + iz * self.l
                    p_ = [x, y, z]
                    new_ltc = Lattice(p_, self.l, n_trd=self.n_trd)
                    R.append(new_ltc)
                P.append(R)
            self.lattrix.append(P)

    def intersect_trident(self, tr:Trident):
        pass

    def intersect_ray(self, ray):
        """
        `ray`: must be a instance of `pyHIFU.ray.Ray` (PowRay / AuxRay)
        return in and out points
        """
        p = list()
        for face in self:
            p_ = ray.intersect_plane(face)
            if (p_ is not None) and ray.has_point(p_):
                p.append(p_)
        return sorted(p, key=lambda xyz: ray.xyz_to_t(xyz))
    
    def affected_lattices(self, ray):
        """ find all affected lattices of ray """
        p_list = self.intersect_ray(ray)
        if len(p_list) == 1:
            # TODO what if ray start inside and only one intersection point
            # extend to list of ray.p
            p_list = [ray.p].extend(p_list)
        assert len(p_list) == 2
        concat_p_list = np.concatenate((p_list[0], p_list[1])).reshape(2,3)
        p_min = concat_p_list.min(axis=0)
        p_max = concat_p_list.max(axis=0)
        start = np.floor((p_min - self.o1) / self.l)
        end = np.ceil((p_max - self.o1) / self.l)
        s = start.astype(int)
        e = end.astype(int)
        return s, e

    def update_all_lattices(self, ray):
        pass
        

class Lattice():
    def __init__(self, p, l, n_trd=256):
        self.o = np.array(p)
        self.a = np.array([l, 0, 0])
        self.b = np.array([0, l, 0])
        self.c = np.array([0, 0, 1])
        # super().__init__(o=o, a=a, b=b, c=c)
        
        self.warehouse = np.zeros((n_trd,))  # to store all the contribution from transducers

    def intersect_ray(self, r):
        pass

if __name__ == "__main__":
    from pyHIFU.geometric.lines import Ray
    r = Ray([0,0,0], [0.1, 0.05, 0.04])
    # r2 = Ray([0.12,0,0], [0.15, 0.05, 0.04])
    B = Box(0.1,-0.1,-0.1, 0.3, 0.1,0.1, 0.01)
    print(B.affected_lattices(r))
    # print(B.affected_lattices(r2))
