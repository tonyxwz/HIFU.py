from pyHIFU.geometric.volumes import *
import unittest
import numpy as np


class Test_CylinderMethods(unittest.TestCase):
    def test_adjmtx(self):
        C = Cylinder([0,0,0], 5, [0,0,2])
        self.assertTrue(np.all(C.adj_mtx == np.array([[1,0,1],
                                                  [0,1,1],
                                                  [1,1,1]])))