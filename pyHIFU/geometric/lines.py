import numpy as np
from pyHIFU.geometric.vec3 import Vec3


class Line(object):
    def __init__(self, p, d):
        self.p = np.array(p)
        self.d = Vec3.normalize(d)
        self.unit_vector = self.d
        self._d_as_assigned = d
        # self.a
        # self.b
        # self.c
        # self.d
        # self.abcd = 

    def __eq__(self, other):
        return (Vec3.are_parallel(self.d, other.d) and
            self.has_point(other.p))

    def __str__(self):
        return str(self.__dict__)

    def to_coordinate(self, t):
        """ convert position on a line to coordinate """
        return self.p + t * self.d

    @property
    def parameters(self):
        return {"p": self.p, "d": self.d}

    def intersect_line(self, line2):
        pass
    
    def intersect_plane(self, plane):
        return plane.intersect_line(line=self)

    def has_point(self, point):
        return all(np.cross(point-self.p, self.d) == 0)

    def is_perpendicular(self, vector=None, line=None):
        if vector is None:
            vector = line.d
        return self.d.dot(vector) == 0

    def is_parallel(self, other):
        return Vec3.are_equal(self.d, other.d)
    
    def find_foot(self, point):
        """ food in perpendicular """
        if type(point) is not np.ndarray:
            point = np.array(point)
        
        t = self.d.dot(point - self.p) / self.d.dot(self.d)
        return self.to_coordinate(t)
        
    def distance_to_point(self, point):
        foot = self.find_foot(point)
        return np.linalg.norm(foot - point)


class Ray(Line):
    def __init__(self, p, d, **kwargs):
        super().__init__(p, d)
        self.start = self.p

    def __eq__(self, other):
        return (Vec3.are_equal(self.d, other.d) and
            Vec3.are_equal(self.start, other.start))
    
    def has_point(self, point):
        if (super().has_point(point) and
            self.d.dot(point - self.start) >= 0):
            return True
        else:
            return False


class Segment(Line):
    def __init__(self, p, d, l=None, **kwargs):
        super().__init__(p, d)
        self.start = self.p
        if l is None:
            self.length = np.linalg.norm(d)
            self.end = self.p + d
        else:
            self.length = l
            self.end = self.to_coordinate(l)

    def __eq__(self, other):
        return (super().__eq__(other) and 
            Vec3.are_equal(self.start, other.start) and
            Vec3.are_equal(self.end, other.end))
        
    def has_point(self, point):
        if (super().has_point(point) and 
            self.d.dot(point - self.start) >= 0 and
            self.d.dot(point - self.end) <= 0):
            return True
        else:
            return False


# Point is just np.ndarray shape = (3,)
# class Point(np.ndarray):
#     def __init__(self, array):
#         super().__init__(array)
#     def on_line(self, line):
#         pass
#     def on_plane(self, plane):
#         pass
