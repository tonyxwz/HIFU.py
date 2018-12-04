import numpy as np
import numpy.linalg as linalg
from sympy.geometry.line import LinearEntity3D, Ray3D, Line3D
from sympy.geometry.line import Ray as SympyRay
from sympy.geometry.line import Line as SympyLine
from sympy.geometry.plane import Plane as SympyPlane


class Vector(object):
    def __init__(self, start, direction, length):
        self.start = start
        self.direction = direction / np.linalg.norm(direction)
        self.end = start + direction * length

    @property
    def length(self):
        return np.linalg.norm(self.end - self.start)


class Plane(SympyPlane):
    """
    Plane class <- sympy.Plane
    """
    @classmethod
    def from_two_vec(self, v1, v2):
        """
        construct plane from two vectors
        """
        self.p = v1.start
        mtx = np.ones([4,4])
        mtx[0:3, 0] = v1.start
        mtx[0:3, 1] = v1.end
        mtx[0:3, 2] = v2.start
        mtx[0:3, 3] = v2.end

        if np.linalg.det(mtx) == 0:
            tmp = np.cross(v1, v2)
            self.n = tmp / np.linalg.norm(tmp)
        else:
            raise Exception("Input vectors are not on the same plane")

    @property
    def get_equation(self):
        """
        get equation of ax + by + cz = d
        """
        abcd = dict()
        abcd["a"] = self.normal_vector[0]
        abcd["b"] = self.normal_vector[1]
        abcd["c"] = self.normal_vector[2]
        abcd["d"] = np.sum(np.multiply(self.p1, self.normal_vector))
        return abcd


    def interp(self, line):
        """
        Intersection point with a line
        """
        if type(line) is not Line:
            line = Line(line.start, line.direction)
        
        A = np.ndarray([3, 3])
        A[:, 1] = self.n
        D = np.ndarray([3,1])
        return A**(-1) * D


    def interl(self, plane):
        """
        Intersection line with a plane
        """
        pass


class PetitePlane(object):
    """ faster and more simplified"""
    def __init__(self, p, n):
        self.p = np.array(p)
        self.normal_vector = np.array(n)

    def intersect_line(self, line):
        """ intersection point with a line"""
        b = self.normal_vector.dot(line.d)
        if b == 0:
            return None
        else:
            t = self.normal_vector.dot(self.p - line.p)
            return line.to_coordinate(t)
    
    def intersect_plane(self, plane2):
        return "line"
    
    def get_matrix(self):
        pass


class Line(object):
    def __init__(self, p, d):
        self.p = np.array(p)
        self.d = np.array(d)

    def to_coordinate(self, t):
        return self.p + t * self.d

    @property
    def parameters(self):
        return {"p": self.p, "d": self.d}

    def intersect_line(self, line2):
        pass
    
    def intersect_plane(self, plane):
        return plane.intersect_line(self)
        

# Point is np.ndarray
# class Point(np.ndarray):
#     def __init__(self, x, y, z):
#         self.x = x
#         self.y = y
#         self.z = z
