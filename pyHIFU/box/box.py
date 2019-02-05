from pyHIFU.geometric.volumes import Cuboid
from pyHIFU.ray import Trident
import numpy as np
from scipy.sparse import coo_matrix
from pyHIFU.geometric.vec3 import EPS, Vec3

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

        self.lx = l # TODO lattice of any width, length, depth
        self.ly = l
        self.lz = l
        self.l = [self.lx, self.ly, self.lz]

        self.n_trd = n_trd  # num of transducers
        self.nx = np.int(np.ceil((xmax-xmin) / self.lx))
        self.ny = np.int(np.ceil((ymax-ymin) / self.ly))
        self.nz = np.int(np.ceil((zmax-zmin) / self.lz))
        self.nxyz = np.array([self.nx, self.ny, self.nz])
        # self.lattrix = list()  # LAT(tice ma)TRIX 3D
        # for ix in range(self.nx):
        #     x = xmin + ix * self.lx
        #     P = list()  # 2D (Plane)
        #     for iy in range(self.ny):
        #         y = ymin + iy * self.ly
        #         R = list()  # 1D (Row)
        #         for iz in range(self.nz):
        #             z = zmin + iz * self.lz
        #             p_ = [x, y, z]
        #             new_ltc = Lattice(p_, self.l, n_trd=self.n_trd)
        #             R.append(new_ltc)
        #         P.append(R)
        #     self.lattrix.append(P)
    
    def intersect_trident(self, tr:Trident):
        pass

    def traversal_along(self, ray, func=None, debug=False):
        """ A implementation of 3DDA algorithm as described in
        "A fast voxel traversal algorithm for ray tracing" 
        by John Amanatides, Andrew Woo, 1984

        http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.42.3443&rep=rep1&type=pdf

        `ray`: pyHIFU.ray.Ray
        """
        st, et = self.rayBoxIntersection(ray)
        if st is not None:
            # 1. determin starting coordinate of the ray: X, Y, Z
            s = ray.to_coordinate(st) - self.o1
            e = ray.to_coordinate(et) - self.o1
            
            s_indices = s // self.l
            e_indices = e // self.l

            ## Compare with the dimensions of the box
            s_indices = np.min([s_indices, self.nxyz - 1], axis=0)
            e_indices = np.min([e_indices, self.nxyz - 1], axis=0)
            ## Compare with 0 (can happen when ray start point cannot be precisely represented)
            s_indices = np.max([s_indices, [0,0,0]], axis=0)
            e_indices = np.max([e_indices, [0,0,0]], axis=0)

            X = int(s_indices[0])
            Y = int(s_indices[1])
            Z = int(s_indices[2])

            X_out = int(e_indices[0])
            Y_out = int(e_indices[1])
            Z_out = int(e_indices[2])

            # 2. stepX, stepY, stepZ
            stepX = 1 if ray.d[0] >= 0 else -1
            stepY = 1 if ray.d[1] >= 0 else -1
            stepZ = 1 if ray.d[2] >= 0 else -1

            # 3. tMaxX, tMaxY, tMaxZ, initial distance from the three directed planes
            ## handle de EPS problem
            m = propermod(s, self.l)

            tMaxX = self.lx - m[0] if stepX == 1 else m[0]
            tMaxY = self.ly - m[1] if stepY == 1 else m[1]
            tMaxZ = self.lz - m[2] if stepZ == 1 else m[2]

            # 4. tDeltaX, tDeltaY, tDeltaZ distance to move on the ray to equal
            # the width/height/depth of a voxel
            # TODO absolute value ?
            tDeltaX = abs(self.lx * ray.d_inv[0])
            tDeltaY = abs(self.ly * ray.d_inv[1])
            tDeltaZ = abs(self.lz * ray.d_inv[2])

            # one dimension goes out of boundary
            while not (X == X_out+stepX or Y == Y_out+stepY or Z == Z_out+stepZ):
                if func:
                    func(X,Y,Z)
                if tMaxX < tMaxY:
                    if tMaxX < tMaxZ:
                        tMaxX = tMaxX + tDeltaX
                        X = X + stepX
                    else:
                        tMaxZ = tMaxZ + tDeltaZ
                        Z = Z + stepZ
                else:
                    if tMaxY < tMaxZ:
                        tMaxY = tMaxY + tDeltaY
                        Y = Y + stepY
                    else:
                        tMaxZ = tMaxZ + tDeltaZ
                        Z = Z + stepZ
            

    def rayBoxIntersection(self, ray):
        """
        find start (tmin) and end (tmax) coordinate of the ray in the lattrix system
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
   

class Lattice():
    def __init__(self, p, l, n_trd=256):
        self.o = np.array(p)
        self.a = np.array([l, 0, 0])
        self.b = np.array([0, l, 0])
        self.c = np.array([0, 0, 1])
        self.center = self.o + (self.a+self.b+self.c) / 2  # for calculating intensity
        self.o2 = self.o + self.a + self.b + self.c

        # super().__init__(o=o, a=a, b=b, c=c)
        
        # self.I_mtx = dict()  # store the total intensity from different bundles
        # self.phi_mtx = dict()  # store total the phase of different bundles
        # self.counter_mtx = dict()  # store how many tridents from different bundles

    def get_initial_values(self, r):
        pass

def propermod(a, b):
    """ return a % b, along with handling EPS problem """
    m = a % b
    for i in range(len(m)):
        if m[i] < EPS or b[i] - m[i] < EPS:
            m[i] = 0
    return m

if __name__ == "__main__":
    from pyHIFU.ray import AuxRay
    r = AuxRay([1.3, 0.57, 0.34], [-0.13, -0.057, -0.034])
    r.end = [0, 0, 0]

    r1 = AuxRay([0, 0, 0], [0.13, 0.057, 0.034])
    r1.end = [1.3, 0.57, 0.34]
    # r2 = Ray([0.12,0,0], [0.15, 0.05, 0.04])
    B = Box(0.1,-0.1,-0.1, 0.3, 0.1,0.1, 0.005)

    print("r")
    B.traversal_along(r, func=print, debug=True)
    print("r1")
    B.traversal_along(r1, func=print, debug=True)
