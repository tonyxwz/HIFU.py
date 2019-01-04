import unittest
import numpy as np
from pyHIFU.geometric.surfaces import Plane, BarrelShell
from pyHIFU.geometric.lines import Line

class Test_PlaneMethods(unittest.TestCase):
    def test_equal(self):
        p1 = Plane([7/3, 7/3, 7/3], [1,1,1])
        p2 = Plane([2,2,3], [2,2,2])
        self.assertTrue(p1 == p2)
        p3 = Plane([8/3, 8/3, 8/3], [1,1,1])
        p4 = Plane([2,3,3], [2,2,2])
        self.assertTrue(p3 == p4)

class Test_BarrelShell(unittest.TestCase):
    def test_intersectline(self):
        B = BarrelShell([0,0,0], [0,0,20], 7, radius2=3)
        l = Line([0,-20, 10],[0,1,0])
        self.assertTrue(np.all(B.intersect_line(l) == [0,-7,10]))
