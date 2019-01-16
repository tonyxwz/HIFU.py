import numpy as np
from .vec3 import Vec3
from .vec3 import EPS

class Line(object):
    def __init__(self, p, d=None, p2=None):
        self.p = np.array(p)
        if d is None:
            d = p2 - self.p
        else:
            p2 = self.p + d
        self.p2 = np.array(p2)
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
        # check can intersect
        vn = np.cross(self.d, line2.d)
        # self and line2 are not parallel
        if np.any(vn):
            v1 = self.p - line2.p
            v2 = self.p2 - line2.p2
            # self and line2 are on the same plane (can intersect)
            if np.dot(vn, v1) == 0 and np.dot(vn, v2) == 0:
                # find the right coor component to use
                d1 = self.d
                d2 = line2.d
                for i in range(len(d1)):
                    if not d1[i] == 0:
                        for j in range(len(d2)):
                            if not d2[j] == 0:
                                if not i == j:
                                    break
                k1 = line2.p[i]*line2.d[j] - line2.p[j]*line2.d[i]
                k2 = self.p[i]*line2.d[j] - self.p[j]*line2.d[i]
                k3 = self.d[i]*line2.d[j] - self.d[j]*line2.d[i]
                t = (k1 - k2 ) / k3
                return self.to_coordinate(t)
        return None

    def intersect_plane(self, plane):
        return plane.intersect_line(line=self)

    def has_point(self, point):
        if point is not None:
            return all(np.cross(point-self.p, self.d) <= EPS)
        else:
            return False

    def is_perpendicular(self, vector=None, line=None):
        if vector is None:
            vector = line.d
        return self.d.dot(vector) == 0

    def is_parallel(self, other):
        return not np.any(np.cross(self.d, other.d))

    def find_foot(self, point):
        """ food in perpendicular """
        if type(point) is not np.ndarray:
            point = np.array(point)

        t = self.d.dot(point - self.p) / self.d.dot(self.d)
        return self.to_coordinate(t)

    def distance_to_point(self, point):
        foot = self.find_foot(point)
        return np.linalg.norm(foot - point)

    def perpendicularDirection(self):
        r_point = np.array([0,0,0])
        for c1 in range(len(self.d)):
            if not self.d[c1] == 0:
                # component 1 is not zero
                for c2 in range(c1+1, len(self.d)):
                    if not self.d[c2] == 0:
                        # component 2 is not zero
                        r_point[c1] = -self.d[c2]
                        r_point[c2] = -self.d[c1]
                        return r_point
                # 2 components are zero, then set "next" component to 1
                r_point[(c1 + 1) % 3] = 1
                return r_point
        return None # sth wrong

    def rotate(self, axis, theta):
        return np.dot(Vec3.rotation_matrix(axis, theta), self)

    



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

    def intersect_plane(self, plane):
        p = super().intersect_plane(plane)
        if p is not None:
            if self.has_point(p):
                return p
        return None


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
        return ((Vec3.are_equal(self.start, other.start) and
                 Vec3.are_equal(self.end, other.end)) or
                (Vec3.are_equal(self.start, other.end) and
                 Vec3.are_equal(self.end, other.start)))

    def has_point(self, point):
        if (super().has_point(point) and
            self.d.dot(point - self.start) >= 0 and
                self.d.dot(point - self.end) <= 0):
            return True
        else:
            return False
