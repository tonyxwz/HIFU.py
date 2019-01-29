from pyHIFU.geometric.volumes import Cuboid
import numpy as np


class Box(Cuboid):
    def __init__(self,
                 xmin, ymin, zmin,
                 xmax, ymax, zmax,
                 l, n_trd=256):
        """
        `[xmin, ymin, zmin]` and `[xmax, ymax, zmax]` defines the diagonal of the Box

        `l`: the edge length of lattice
        """
        o = [xmin, ymin, zmin]
        a = [xmax-xmin, 0, 0]
        b = [0, ymax-ymin, 0]
        c = [0, 0, zmax-zmin]
        super().__init__(o, a, b, c)  # o1, o2

        self.l = l
        self.n_trd = 256
        nx = abs(xmax-xmin) // self.l
        ny = abs(ymax-ymin) // self.l
        nz = abs(zmax-zmin) // self.l
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
            
                
                    


class Lattice(Cuboid):
    def __init__(self, p, l, n_trd=256):
        o = np.array(p)
        a = np.array([l, 0, 0])
        b = np.array([0, l, 0])
        c = np.array([0, 0, 1])
        super().__init__(o=o, a=a, b=b, c=c)
        
        self.warehouse = np.zeros((n_trd,))  # to store all the contribution fron transducers
