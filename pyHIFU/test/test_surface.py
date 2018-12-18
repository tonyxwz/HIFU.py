import unittest
import numpy as np
from pyHIFU.geometric.surfaces import Plane

class Test_PlaneMethods(unittest.TestCase):
    def test_equal(self):
        p1 = Plane([7/3, 7/3, 7/3], [1,1,1])
        p2 = Plane([2,2,3], [2,2,2])
        self.assertTrue(p1 == p2)
        p3 = Plane([8/3, 8/3, 8/3], [1,1,1])
        p4 = Plane([2,3,3], [2,2,2])
        self.assertTrue(p3 == p4)