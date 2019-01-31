import numpy as np
from .vec3 import Vec3
from .lines import Line, Ray, Segment
from .curves import Ring


class Plane(object):
    """ 
    All edges should be counter-clockwise
    right hand, cross product direction
    """

    def __init__(self, p, n, **kw):
        self.p = np.array(p)
        self.normal_vector = Vec3.normalize(n)
        self.n = self.normal_vector
        self._p_as_assigned = p
        self._n_as_assigned = n
        self.abcd = np.ones(4)
        self.abcd[0:3] = n
        self.abcd[3] = -self.p.dot(n)
        self.edges = []

    def __eq__(self, other):
        return (Vec3.are_perpendicular(self.p - other.p, self.normal_vector) and
                Vec3.are_equal(self.normal_vector, other.normal_vector))

    def intersect_line(self, line=None, p0=None, vd=None, require_t=False):
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
            if require_t:
                return t
            else:
                return line.to_coordinate(t)

    def intersect_plane(self, plane2):
        return "line"

    def has_point(self, point=None):
        if point is not None:
            point = np.array(point)
            return Vec3.are_perpendicular(self.normal_vector, point - self.p)
        else:
            return False

    def has_line(self, line=None):
        return (self.has_point(point=line.p) and
                line.is_perpendicular(vector=self.normal_vector))

    def n_common_edge(self, other):
        """ return number of common edges """
        ans = 0
        for e1 in self.edges:
            for e2 in other.edges:
                ans += int(e1 == e2)
        return ans

    def distance_to_point(self, point):
        v1 = self.p - point
        d = np.abs(np.dot(v1, self.normal_vector))
        return d


class Circle(Plane):
    def __init__(self, center, radius, normal_vector, radius2=None, **kw):
        # TODO: hollowed circle
        super().__init__(center, normal_vector)
        self.center = self.p
        self.radius = radius
        self.edges = [Ring(self.center, self.radius, self.normal_vector)]
        if radius2 is not None:
            self.radius2 = radius2
            self.edges.append(Ring(self.center, self.radius2, self.normal_vector))
        else:
            self.radius2 = 0
        if "index" in kw:
            self.index = kw['index']
        
    @property
    def area(self):
        return np.pi * np.power(self.radius, 2)

    @property
    def circumference(self):
        return 2 * np.pi * self.radius

    def __eq__(self, other):
        return (all(self.center == other.center) and
                all(self.normal_vector == other.normal_vector) and
                self.radius == other.radius and
                self.radius2 == other.radius2)

    def has_point(self, point):
        # a point is in the circle if the point is on th plane and
        # the distance to origin is less than radius
        # TODO: support for ellipse
        return (super().has_point(point) and
                np.linalg.norm(np.array(point)-self.center) <= self.radius and
                np.linalg.norm(np.array(point)-self.center) >= self.radius2)


class Polygon(Plane):
    def __init__(self, p, vertices, edges, **kw):
        self.vertices = vertices
        self.edges = edges
    
    def __is_polygon(self):
        """ prove the args of `__init__` can form a polygon """
        pass

    def has_point(self, p):
        pass

    def has_segment(self, s):
        pass
    
    def is_containing(self, other):
        """
        if all the vertices of `other` are inside `self` AND
        all the edges has no intersection point, then self is containing other.
        """
        ans = True
        for v in other.vertices:
            ans = ans and self.has_point(v)
        for e in other.edges:
            ans = ans and self.has_line(e)
        return ans


class Triangle(Polygon):
    """ most simple polygon """
    pass


class Rectangle(Plane):
    def __init__(self, p, va, vb, **kw):
        super().__init__(p, np.cross(va, vb))
        self.va = np.array(va)
        self.vb = np.array(vb)
        if "index" in kw:
            self.index = kw['index']

        self.edges.append(Segment(self.p, self.va))
        self.edges.append(Segment(self.p+self.va, self.vb))
        self.edges.append(Segment(self.p+self.va+self.vb, np.negative(self.va)))
        self.edges.append(Segment(self.p+self.vb, np.negative(self.vb)))

        self.vertices = [self.p, self.p+self.va,
                         self.p+self.va+self.vb, self.p+self.vb]

    def __str__(self):
        return str(self.__dict__)

    def __eq__(self, other):
        """ equality regardless of normal vectors """
        ans = True
        for e1 in self.edges:
            q = False
            for e2 in other.edges:
                q = q or (e1 == e2)
            ans = ans and q
        return ans
    def is_overlapping(self, other):
        """ if all the vertices of one rectangle are inside another rect"""
        pass

    def has_point(self, point):
        f1 = self.edges[0].find_foot(point)
        f2 = self.edges[1].find_foot(point)
        if self.edges[0].has_point(f1) and self.edges[1].has_point(f2):
            return True
        else:
            return False

    def intersect_line(self, line=None, p0=None, vd=None, require_t=False):
        p = super().intersect_line(line=line, p0=p0, vd=vd, require_t=require_t)
        if p is not None:
            if self.has_point(p):
                return p
        return None

    def has_rect(self, other):
        pass


class Sphere(object):
    """ Sphere with direction and angle """

    def __init__(self, center, radius, v_axis, angle=np.pi):
        self.center = np.array(center)
        self.radius = np.float(radius)
        if angle > 0 and angle <= np.pi:
            self.angle = angle
        else:
            raise ValueError()
        self.v_axis = np.array(v_axis)

    def has_point(self, point):
        v = point - self.center
        dist = np.linalg.norm(v)
        if (dist - self.radius) < np.finfo(float).eps:
            costheta = np.dot(self.v_axis, v) / (dist * self.radius)
            if costheta < np.cos(self.angle):
                return True
        return False

    def tan_plane_at(self, point):
        vn = point - self.center
        return Plane(point, vn)

    def r_to_xyz(self, r):
        # intersection point of vector r with self
        return


class BarrelShell(object):
    """ Cylinder surface """

    def __init__(self, center, v_axis, radius, **kw):
        self.v_axis = np.array(v_axis)
        self.radius = radius
        # here center is on the bottom because this cylinder is directed
        self.center = center
        self.length = np.linalg.norm(self.v_axis)
        self.unit_vector = Vec3.normalize(self.v_axis)
        self.axis = Segment(center, v_axis)
        self.edges = [Ring(self.center, self.radius, self.unit_vector),
                      Ring(self.center + self.v_axis, self.radius, self.unit_vector)]
        if "index" in kw.keys():
            self.index = kw["index"]

    def has_point(self, point):
        if type(point) is not np.ndarray:
            point = np.array(point)
        foot = self.axis.find_foot(point)
        dist = self.axis.distance_to_point(point)
        return (dist - self.radius < np.finfo(float).eps and
                self.axis.has_point(foot))

    def tan_plane_at(self, point):
        # TODO: calculate tangential plane for acoustics calculation
        l = Line(self.center, self.axis)
        foot = l.find_foot(point)
        vn = point - foot
        return Plane(point, vn)
        
    def n_common_edge(self, other):
        """ return number of common edges """
        ans = 0
        for e1 in self.edges:
            for e2 in other.edges:
                ans += int(e1 == e2)
        return ans

    def __eq__(self, other):
        if type(other) is BarrelShell:
            return (self.radius == other.radius and 
                    self.axis == other.axis)
        else: return False

    def intersect_line(self, line):
        vn = np.cross(line.d, self.unit_vector)
        pl = Plane(line.p, vn)
        d = pl.distance_to_point(self.center)
        # l2: 
        l2 = Segment(self.center+pl.normal_vector*d, self.v_axis)
        if not pl.has_line(l2):
            l2 = Segment(self.center-pl.normal_vector*d, self.v_axis)

        p_center = l2.intersect_line(line)
        d2 = np.sqrt(self.radius**2 - d**2)
        costheta = np.dot(self.unit_vector, line.d)
        sintheta = np.sqrt(1-costheta**2)
        t = np.abs(d2 / sintheta)
        p1 = p_center - line.unit_vector * t
        return p1

