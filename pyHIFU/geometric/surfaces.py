import numpy as np
from pyHIFU.geometric.vec3 import Vec3
from pyHIFU.geometric.lines import Line, Ray, Segment


class Plane(object):
    """ faster and more simplified """

    def __init__(self, p, n):
        self.p = np.array(p)
        self.normal_vector = Vec3.normalize(n)
        self._p_as_assigned = p
        self._n_as_assigned = n
        self.abcd = np.ones(4)
        self.abcd[0:3] = n
        self.abcd[3] = -self.p.dot(n)

    def __eq__(self, other):
        return (Vec3.are_perpendicular(self.p - other.p, self.normal_vector) and
                Vec3.are_equal(self.normal_vector, other.normal_vector))

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

    def has_point(self, point=None):
        point = np.array(point)
        return Vec3.are_perpendicular(self.normal_vector, point - self.p)

    def has_line(self, line=None):
        return (self.has_point(point=line.p) and
                line.is_perpendicular(vector=self.normal_vector))

    def get_matrix(self):
        pass


class Circle(Plane):
    def __init__(self, center, radius, normal_vector, radius2=None):
        super().__init__(center, normal_vector)
        self.center = self.p
        self.radius = radius
        if radius2 is not None:
            self.radius2 = radius2
        else:
            self.radius2 = radius

    @property
    def area(self):
        return np.pi * np.power(self.radius, 2)

    @property
    def circumference(self):
        return 2 * np.pi * self.radius

    @property
    def edge(self):
        # TODO: return the a shape of circle? really needed?
        pass

    def __eq__(self, other):
        return (all(self.center == other.center) and
                all(self.normal_vector == other.normal_vector) and
                self.radius == other.radius)

    def has_point(self, point):
        # a point is in the circle if the point is on th plane and
        # the distance to origin is less than radius
        # TODO: support for ellipse
        return (super().has_point(point) and
                np.linalg.norm(np.array(point)-self.center) < self.radius)


class Triangle(Plane):
    pass


class Polygon(Plane):
    pass


class Rectangle(Plane):
    def __init__(self, p, va, vb):
        super().__init__(p, np.cross(va, vb))
        self.va = np.array(va)
        self.vb = np.array(vb)

    def has_point(self, point):
        diag = np.abs(self.va + self.vb)
        if super().has_point(point) and all(diag >= np.abs(point-self.p)):
            return True
        else:
            return False


class Sphere(object):
    """ Sphere with direction and angle """
    def __init__(self, center, radius, axis_vector, angle=np.pi):
        self.center = np.array(center)
        self.radius = np.float(radius)
        if angle > 0 and angle <= np.pi:
            self.angle = angle
        else:
            raise ValueError()
        self.axis_vector = np.array(axis_vector)

    def has_point(self, point):
        v = point - self.center
        dist = np.linalg.norm(v)
        if (dist - self.radius) < np.finfo(float).eps:
            costheta = np.dot(self.axis_vector, v) / (dist * self.radius)
            if costheta < self.angle:
                return True
        return False


class Barrel(object):
    """ Cylinder surface """

    def __init__(self, axis_vector, radius, center):
        self.axis_vector = np.array(axis_vector)
        self.radius = radius
        # here center is on the bottom because this cylinder is directed
        self.center = center
        self.length = np.linalg.norm(self.axis_vector)
        self.unit_vector = Vec3.normalize(self.axis_vector)
        self.axis = Segment(center, axis_vector)

    def has_point(self, point):
        if type(point) is not np.ndarray:
            point = np.array(point)
        foot = self.axis.find_foot(point)
        dist = self.axis.distance_to_point(point)
        return (dist - self.radius < np.finfo(float).eps and
                self.axis.has_point(foot))
