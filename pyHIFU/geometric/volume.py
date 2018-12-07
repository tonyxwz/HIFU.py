import numpy as np
from pyHIFU.geometric.vec3 import Vec3
from pyHIFU.geometric.lines import Line, Ray, Segment
from pyHIFU.geometric.surfaces import Plane, Circle, Sphere, Rectangle
# a volume is a combination of different surfaces.


class Volume(list):
    def __init__(self, face_list, adj_mtx=None, config_file=None):
        pass


class Ball(Volume):
    def __init__(self, radius=None, normal_vector=None, angle=np.pi):
        pass


class Cylinder(Volume):
    def __init__(self, radius=None, axis_vector=None, length=None):
        pass


class Cuboid(Volume):
    def __init__(self, o, a, b, c):
        self.o1 = np.array(o)
        self.a = np.array(a)
        self.b = np.array(b)
        self.c = np.array(c)
        self.o2 = self.o1 + self.a + self.b + self.c


        self.face_ab_1 = Rectangle(self.o1, self.a, self.b)
        self.face_bc_1 = Rectangle(self.o1, self.b, self.c)
        self.face_ac_1 = Rectangle(self.o1, self.a, self.c)
        self.face_ab_2 = Rectangle(self.o2, -self.a, -self.b)
        self.face_bc_2 = Rectangle(self.o2, -self.b, -self.c)
        self.face_ac_2 = Rectangle(self.o2, -self.a, -self.c)

        face_list = [self.face_ab_1, self.face_bc_1, self.face_ac_1,
                     self.face_ab_2, self.face_bc_2, self.face_ac_2]

        adj_mtx = []
