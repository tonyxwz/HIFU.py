import numpy as np
from numba import vectorize
import unittest
from pyHIFU.transducer import TElement, Transducer
from pyHIFU.geometric.vec3 import EPS


@vectorize(['float32(float32, float32)'], target='cuda')
def Add(a, b):
    return a + b

class Test_NumbaCuda(unittest.TestCase):
    def test_add(self):
        # Initialize arrays
        N = 100000
        A = np.ones(N, dtype=np.float32)
        B = np.ones(A.shape, dtype=A.dtype)
        C = np.empty_like(A, dtype=A.dtype)
        # Add arrays on GPU
        C = Add(A, B)
        self.assertTrue(np.all(C == A*2))

class Test_TransducerElement(unittest.TestCase):
    

    def test_buildtridents(self):
        tcenter = [-14, 0, 0]
        elid = 0
        te = TElement(tcenter, elid)
        for tr in te:
            a =np.dot(tr.pow_ray.d, te.axial_ray.d)
            costheta = a / (np.linalg.norm(tr.pow_ray.d) * np.linalg.norm(te.axial_ray.d))
            self.assertTrue(costheta >= np.sqrt(2)/2)

    def test_tridentangle(self):
        tcenter = [-14, 0, 0]
        elid = 0
        te = TElement(tcenter, elid)
        for tr in te:
            costheta = np.dot(tr.pow_ray.d, tr.aux_ray1.d) / (np.linalg.norm(tr.pow_ray.d) *
                                                             np.linalg.norm(tr.aux_ray1.d))
            self.assertTrue(costheta - np.cos(1e-4) < EPS)
            n1 = np.cross(tr.pow_ray.d, tr.aux_ray1.d)
            n2 = np.cross(tr.pow_ray.d, tr.aux_ray2.d)
            self.assertTrue(np.abs(np.dot(n1, n2)) < EPS)

