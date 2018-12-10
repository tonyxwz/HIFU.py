import numpy as np
from pyHIFU.geometric.vec3 import Vec3
from pyHIFU.geometric.lines import Line, Ray, Segment
from pyHIFU.geometric.surfaces import Plane, Circle, Sphere, Rectangle
# a volume is a combination of different surfaces.


class Volume(list):
    def __init__(self, sides_list, adj_mtx=None, config_file=None):
        super().__init__(sides_list)
        self.adj_mtx = adj_mtx

    def common_edge(self, id1, id2):
        return self.adj_mtx[id1, id2]


class Ball(Volume):
    def __init__(self, radius=None, normal_vector=None, angle=np.pi):
        pass


class Cylinder(Volume):
    def __init__(self, radius, v_axis, length, radius2=None, v_axis2=None):
        pass


class Cube(Volume):
    def __init__(self, o, a, b, c):
        self.o1 = np.array(o)
        self.a = np.array(a)
        self.b = np.array(b)
        self.c = np.array(c)
        self.o2 = self.o1 + self.a + self.b + self.c

        sides_list = [
            Rectangle(self.o1, self.a, self.c),
            Rectangle(self.o1, self.b, self.a),
            Rectangle(self.o1, self.c, self.b),
            Rectangle(self.o2, np.negative(self.a), np.negative(self.b)),
            Rectangle(self.o2, np.negative(self.b), np.negative(self.c)),
            Rectangle(self.o2, np.negative(self.c), np.negative(self.a))]

        adj_mtx = np.eye(6)
        for i, s1 in enumerate(sides_list):
            for j, s2 in enumerate(sides_list):
                if s1.n_common_edge(s2):
                    adj_mtx[i, j] = 1

        super().__init__(sides_list, adj_mtx=adj_mtx)

    def __str__(self):
        return str([self.o1, self.a, self.b, self.c])


if __name__ == '__main__':
    c = Cube([0, 0, 0], [5, 0, 0], [0, 3, 0], [0, 0, 4])
    print(c)
