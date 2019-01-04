import unittest
import numpy as np
from pyHIFU.geometric.lines import *
from pyHIFU.geometric.surfaces import Plane
from pyHIFU.geometric.vec3 import Vec3


class Test_LineMethods(unittest.TestCase):
    def test_equal(self):
        l1 = Line([1,1,0], p2=[0,2,1])
        l2 = Line([0.5, 1.5, 0.5], d=[1,-1,-1])
        self.assertTrue(l1 == l2)
    
    def test_intersect_plane(self):
        plane = Plane([7/4,0,0], [4,5,3])
        l = Line([1,1,0], d=[-1,1,1])
        self.assertTrue(Vec3.are_equal(l.intersect_plane(plane), np.array([1.5,0.5,-0.5])))

    def test_has_point(self):
        l = Line([6,5,9], d=[-5,-4,-8])
        self.assertTrue(l.has_point([1,1,1]))

    def test_perpendicular(self):
        l = Line([6,5,9], d=[-5,-4,-8])
        l2 = Line([7,7,7], d=[-5,-4,-8])
        self.assertTrue(l.is_parallel(l2))

    def test_find_foot(self):
        l = Line([1,1,1],[3,4,0])
        p = l.find_foot([1,6,1])
        self.assertTrue(np.all(np.equal(p, np.array([12/5 + 1, 16/5 + 1, 1]))))

    def test_distance_to_point(self):
        l = Line([1,1,1],[-1,1,0])
        self.assertEqual(l.distance_to_point([0,0,0]), np.sqrt(3))


class Test_SegmentMethods(unittest.TestCase):
    def test_equal(self):
        seg1 = Segment([2,3,4], [1,1,1])
        seg2 = Segment([3,4,5], [-1,-1,-1])
        self.assertTrue(seg1 == seg2)

    def test_has_point(self):
        seg1 = Segment([2,3,4], [1,1,1])
        self.assertFalse(seg1.has_point([3.1,4.1,5.1]))
        self.assertTrue(seg1.has_point([2.5,3.5,4.5]))


