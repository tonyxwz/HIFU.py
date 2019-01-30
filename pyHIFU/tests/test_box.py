import unittest
from pyHIFU.geometric.lines import Ray
from pyHIFU.box import Box
import numpy as np
from pyHIFU.geometric.vec3 import Vec3

class Test_BoxMethods(unittest.TestCase):
    def test_intersect_ray(self):
        p1 = [2, 4, 6]
        p2 = [5, 4, 6]
        d = [1.5, 1.8, 0.9]
        r = Ray(p1, d)
        r2 = Ray(p2, d)  # two parallel rays
        
        p3 = [3.5, 0, 0]
        p4 = [8, 9, 9]
        B = Box(*p3, *p4, 1)
        plist1 = B.intersect_ray(r)
        plist2 = B.intersect_ray(r2)

        self.assertTrue(len(plist1) == 2)
        self.assertTrue(len(plist2) == 1)
        
        p3 = [3.5, -20, -20]
        p4 = [5, 20, 20]
        B = Box(*p3, *p4, 1)

        p5 = [0, 1, 1]
        d1 = [2, 3, -2]
        d2 = [3, -1.1, 1]
        r3 = Ray(p5, d1)
        r4 = Ray(p5, d2)  # two rays from the same start

        plist3 = B.intersect_ray(r3)
        plist4 = B.intersect_ray(r4)

        ratio1 = (plist4[0] - p5) / (plist4[1] - p5)
        ratio2 = (plist3[0] - p5) / (plist3[1] - p5)
        self.assertTrue(Vec3.are_equal(ratio1, ratio2))

    
