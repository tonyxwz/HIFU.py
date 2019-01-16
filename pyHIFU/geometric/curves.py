import numpy as np
from .vec3 import Vec3


class Ring(object):
    def __init__(self, c, r, n, radius_angle0=None):
        self.c = np.array(c)
        self.r = np.float(r)
        self.n =np.array(n)
        # the radial direction whose reference angle is 0
        if radius_angle0 is not None:
            self.radius_angle0 = np.array(radius_angle0)

    def __eq__(self, other):
        if type(other) is Ring:
            return (all(self.c == other.c) and
                    (self.r == other.r) and
                    Vec3.are_parallel(self.n, other.n))
        return False

    def angletoxyz(self, angle):
        # convert angle to coordinates
        pass

    def getPerpendicularDirection(self):
        r_point = np.array([0,0,0])
        for c1 in range(len(self.n)):
            if not self.n[c1] == 0:
                # component 1 is not zero
                for c2 in range(c1+1, len(self.n)):
                    if not self.n[c2] == 0:
                        # component 2 is not zero
                        r_point[c1] = -self.n[c2]
                        r_point[c2] = -self.n[c1]
                        return r_point
                # 2 components are zero, then set "next" component to 1
                r_point[(c1 + 1) % 3] = 1
                return r_point
        return None # sth wrong

                    
