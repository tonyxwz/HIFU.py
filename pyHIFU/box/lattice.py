import numpy as np


class Lattice():
    def __init__(self, p, l, n_trd=256):
        self.o = np.array(p)
        self.a = np.array([l[0], 0, 0])
        self.b = np.array([0, l[1], 0])
        self.c = np.array([0, 0, l[2]])
        self.center = self.o + (self.a+self.b+self.c) / 2  # for calculating intensity
        self.o2 = self.o + self.a + self.b + self.c

        # super().__init__(o=o, a=a, b=b, c=c)
        
        # self.I_mtx = dict()  # store the total intensity from different bundles
        # self.phi_mtx = dict()  # store total the phase of different bundles
        # self.counter_mtx = dict()  # store how many tridents from different bundles

    def get_initial_values(self, r):
        pass

