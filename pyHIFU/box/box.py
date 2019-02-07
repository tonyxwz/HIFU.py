import numpy as np
from scipy.sparse import coo_matrix

from pyHIFU.geometric.vec3 import EPS, Vec3
from pyHIFU.geometric.volumes import Cuboid
from pyHIFU.ray import Trident
from pyHIFU.visualization.mkplots import plot_box, plot_lattice, plot_ray
from pyHIFU.box.lattice import Lattice


class Box(Cuboid):
    def __init__(self, x1, y1, z1, x2, y2, z2,
                 l=None, lx=None, ly=None, lz=None, n_trd=256):
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

        self.lx = lx or l
        self.ly = ly or l
        self.lz = lz or l
        self.l = np.array([self.lx, self.ly, self.lz])

        self.n_trd = n_trd  # num of transducers
        self.nx = np.int(np.round((xmax-xmin) / self.lx))
        self.ny = np.int(np.round((ymax-ymin) / self.ly))
        self.nz = np.int(np.round((zmax-zmin) / self.lz))
        self.nxyz = np.array([self.nx, self.ny, self.nz])

        self.abc = self.nxyz * self.l  # the length of edges (after rounding)

        self.lattrix = list()  # LAT(tice ma)TRIX 3D
        for ix in range(self.nx):
            x = xmin + ix * self.lx
            P = list()  # 2D (Plane)
            for iy in range(self.ny):
                y = ymin + iy * self.ly
                R = list()  # 1D (Row)
                for iz in range(self.nz):
                    z = zmin + iz * self.lz
                    p_ = [x, y, z]
                    new_ltc = Lattice(p_, self.l, n_trd=self.n_trd)
                    R.append(new_ltc)
                P.append(R)
            self.lattrix.append(P)
    
    def intersect_trident(self, tr:Trident):
        pass

    def traversal_along(self, ray, func=None, debug=False, ax=None, color=[]):
        """ A implementation of 3DDA algorithm as described in
        "A fast voxel traversal algorithm for ray tracing" 
        by John Amanatides, Andrew Woo, 1984

        http://www.cse.yorku.ca/~amana/research/grid.pdf

        `ray`: pyHIFU.ray.Ray
        """
        st, et = self.rayBoxIntersection(ray)
        if st is not None:
            # 2. stepX, stepY, stepZ
            stepX = 1 if ray.d[0] >= 0 else -1
            stepY = 1 if ray.d[1] >= 0 else -1
            stepZ = 1 if ray.d[2] >= 0 else -1

            # 3. tDeltaX, tDeltaY, tDeltaZ distance to move on the ray to equal
            # the width/height/depth of a voxel
            tDeltaX = abs(self.lx * ray.d_inv[0])
            tDeltaY = abs(self.ly * ray.d_inv[1])
            tDeltaZ = abs(self.lz * ray.d_inv[2])

            # 1. determin starting coordinate of the ray: X, Y, Z
            s0 = ray.to_coordinate(st)
            e0 = ray.to_coordinate(et)
            if ax:
                ax.scatter([s0[0], e0[0]], [s0[1], e0[1]], [s0[2], e0[2]], color=color)
            
            s = s0 - self.o1  # subtraction -> float risk
            e = e0 - self.o1

            # Quotient, Remainder of s(start) and e(end)
            qs, rs = self.properDivide(s, self.l, [stepX, stepY, stepZ])
            if np.any(qs < 0) or np.any(qs > self.nxyz-1):
                # the ray start at on face of the box and is leaving the box
                return
            # if function didn't return from calculating the start point,
            # there shouldn't be a problem calculation the end point.
            qe, re = self.properDivide(e, self.l, [-stepX, -stepY, -stepZ])

            X = int(qs[0])
            Y = int(qs[1])
            Z = int(qs[2])
            
            X_out = int(qe[0])
            Y_out = int(qe[1])
            Z_out = int(qe[2])

            tMaxX = (1 - rs[0] / self.lx) * tDeltaX if stepX > 0 else rs[0] / self.lx * tDeltaX
            tMaxY = (1 - rs[1] / self.ly) * tDeltaY if stepY > 0 else rs[1] / self.ly * tDeltaY
            tMaxZ = (1 - rs[2] / self.lz) * tDeltaZ if stepZ > 0 else rs[2] / self.lz * tDeltaZ

            # tMaxX = self.lx - m[0] if stepX == 1 else m[0]
            # tMaxY = self.ly - m[1] if stepY == 1 else m[1]
            # tMaxZ = self.lz - m[2] if stepZ == 1 else m[2]

            # one dimension goes out of boundary
            while not (X == X_out+stepX or Y == Y_out+stepY or Z == Z_out+stepZ):
                if func:
                    func(X,Y,Z)
                if ax:
                    if not len(color):
                        color = np.random.random((1,3))
                    plot_lattice(X, Y, Z, self.l, self.o1, ax, color=color)
                if debug:
                    print("tMaxX, tMaxY, tMaxZ:", tMaxX, tMaxY, tMaxZ)
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

    def properDivide(self, a, b, s):
        """ return quotient and remainder, along with handling EPS problem
        `a`: upper
        `b`: lower
        `d`: deltas
        `s`: steps
        """
        r = a % b
        q = a // b
        for i in range(len(r)): # x y z
            if r[i] < EPS and r[i] >= 0:
                if s[i] < 0:
                    q[i] = q[i] - 1
                    r[i] = b[i]
                else:  # step >= 0
                    r[i] = 0
            elif b[i] - r[i] < EPS:
                if s[i] > 0:
                    q[i] = q[i] + 1
                    r[i] = b[i]
                else:  # step >= 0
                    r[i] = 0
        return q, r


if __name__ == "__main__":
    # testing
    from pyHIFU.ray import AuxRay
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn
    from mpl_toolkits.mplot3d import Axes3D
    import mpl_toolkits.mplot3d as mp3d
    from random import random
    import time

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.set_aspect("equal")

    t0 = time.time()
    B = Box(0.1,-0.1,-0.1, 0.3, 0.1,0.1, 0.005)

    r = AuxRay([1.3, 0.57, 0.34], [-0.13, -0.057, -0.034])
    r.end = [0, 0, 0]
    plot_ray(r, ax)
    print("r")
    B.traversal_along(r, func=print, debug=True, ax=ax, color=[1,0,0])

    r1 = AuxRay([0, 0, 0], [0.13, 0.057, 0.034])
    r1.end = [1.3, 0.57, 0.34]
    plot_ray(r1, ax, linestyle="--", color='g')
    print("r1")
    B.traversal_along(r1, func=print, debug=True, ax=ax, color=[0,0,1])

    r2 = AuxRay([0.02, 0.05, 0.01], [0.17, 0.057, 0.034])
    r2.endt = 1
    plot_ray(r2, ax, color='g')
    print("r2")
    B.traversal_along(r2, func=print, debug=True, ax=ax, color=[0,0,1,0.5])

    r3 = AuxRay([0.87 , 0.335, 0.18 ], [-0.17, -0.057, -0.034])
    r3.endt = 1
    plot_ray(r3, ax, linestyle="--", color='r')
    print("r3")
    B.traversal_along(r3, func=print, debug=True, ax=ax, color=[0,1,0,0.5])
    print(time.time() - t0)
    plot_box(B, ax, title="Box intersection")
    plt.show()
