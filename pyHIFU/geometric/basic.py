import numpy as np


class Plane(object):
    """ faster and more simplified """
    def __init__(self, p, n):
        self.p = np.array(p)
        self.normal_vector = np.array(n)

    def intersect_line(self, line=None, p0=None, vd=None):
        """ intersection point with a line """
        if line is not None:
            if (p0 is not None) or (vd is not None):
                raise Exception("Only need to define line once")
        elif (p0 is not None) and (vd is not None):
            line = Line(np.array(p0), np.array(vd))
        else:
            raise Exception("illegal argument")
        
        b = self.normal_vector.dot(line.d)
        if b == 0:
            return None
        else:
            t = self.normal_vector.dot(self.p - line.p) / b
            return line.to_coordinate(t)
    
    def intersect_plane(self, plane2):
        return "line"
    
    def get_matrix(self):
        pass


class Line(object):
    def __init__(self, p, d):
        self.p = np.array(p)
        self.d = np.array(d)
        self.unit_vector = self.d / np.linalg.norm(self.d)

    def to_coordinate(self, t):
        return self.p + t * self.d

    @property
    def parameters(self):
        return {"p": self.p, "d": self.d}

    def intersect_line(self, line2):
        pass
    
    def intersect_plane(self, plane):
        return plane.intersect_line(line=self)
        
class Ray(Line):
    def __init__(self, p, d, **kwargs):
        super().__init__(p, d)



# Point is just np.ndarray shape = (3,)
# class Point(np.ndarray):
#     def __init__(self, x, y, z):
#         self.x = x
#         self.y = y
#         self.z = z
