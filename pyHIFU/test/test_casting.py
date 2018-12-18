import numpy as np
from numba import vectorize
import unittest


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
