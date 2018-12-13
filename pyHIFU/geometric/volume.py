import numpy as np
from .vec3 import Vec3
from .lines import Line, Ray, Segment
from .surfaces import Plane, Circle, Sphere, Rectangle


class Volume(list):
    """
    General class for all geometric volumes such as ball, cylinder, cube
    not patient volume yet
    """
    def __init__(self, sides_list, adj_mtx=None, shape_dict=None):
        super().__init__(sides_list)
        self.adj_mtx = adj_mtx

    def common_edge(self, id1, id2):
        return self.adj_mtx[id1, id2]

    def adj_at(self, volume2):
        """
        return the index of the sides of the face
        where self and volume2 is adjacent
        """
        # TODO: judge two planes, coplanary but not completely the same
        r = list()
        for f1 in self:
            for f2 in volume2:
                if f1 == f2:
                    r.append(f1.index)
        return r
        


class Ball(Volume):
    def __init__(self, radius=None, normal_vector=None, angle=np.pi):
        pass


class Cylinder(Volume):
    def __init__(self, radius, v_axis, length, radius2=None, v_axis2=None):
        pass


class Cuboid(Volume):
    def __init__(self, o=None, a=None, b=None, c=None):
        # print("Cuboid::kw:", kw)
        # print("Cuboid::a,b,c:", a, b, c)
        self.o1 = np.array(o)
        self.a = np.array(a)
        self.b = np.array(b)
        self.c = np.array(c)
        self.o2 = self.o1 + self.a + self.b + self.c

        sides_list = [
            Rectangle(self.o1, self.a, self.c, index=0),
            Rectangle(self.o1, self.b, self.a, index=1),
            Rectangle(self.o1, self.c, self.b, index=2),
            Rectangle(self.o2, np.negative(self.a), np.negative(self.b), index=3),
            Rectangle(self.o2, np.negative(self.b), np.negative(self.c), index=4),
            Rectangle(self.o2, np.negative(self.c), np.negative(self.a), index=5)]
        #adjacency matrix of volume, items indication common edge
        adj_mtx = np.eye(6)
        # TODO calculate only half of the matrx because symmetry property
        for i, s1 in enumerate(sides_list):
            for j, s2 in enumerate(sides_list):
                if s1.n_common_edge(s2):
                    adj_mtx[i, j] = 1

        super().__init__(sides_list, adj_mtx=adj_mtx)

    def __str__(self):
        return str([self.o1, self.a, self.b, self.c])

    def common_faces(self, other):
        # TODO return all the common faces of self and other
        pass


if __name__ == '__main__':
    c = Cuboid([0, 0, 0], [5, 0, 0], [0, 3, 0], [0, 0, 4])
    print(c)
