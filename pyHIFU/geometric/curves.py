import numpy as np
from .vec3 import Vec3


class Ring(object):
    def __init__(self, c, r, n):
        self.c = np.array(c)
        self.r = np.float(r)
        self.n =np.array(n)

    def __eq__(self, other):
        if type(other) is Ring:
            return (all(self.c == other.c) and
                    (self.r == other.r) and
                    Vec3.are_parallel(self.n, other.n))
        return False