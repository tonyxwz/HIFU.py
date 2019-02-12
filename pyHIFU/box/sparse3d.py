"""
Implementing sparse matrix in 3d for pyHIFU project
The only operation needed in sampling box is adding

"""
import numpy as np
import cmath

# TODO Tuples can be used as dictionary key because it is inmutable
class Sparse3D():
    def __init__(self, shape):
        self.shape = shape
        self.__data = {}
        self.__default = 0.0

    def __setitem__(self, index, value):
        """ set value to position given in index, where index is a tuple. """
        if index > (0, 0, 0) and index < self.shape:
            self.__data[index] = float(value)
        else:
            raise Exception("Index out of range of sparse matrix")

    def __getitem__(self, index):
        if index > (0, 0, 0) and index < self.shape:
            return self.__data.get(index, self.__default)
        else:
            raise Exception("Index out of range of sparse matrix")

    def cexp(self):
        out = self.__class__(self.shape)
        for k, v in self.__data.items():
            out[k] = cmath.exp(v*1j)
        return out
    
    def __add__(self, other):
        if self.shape == other.shape:
            out = self.__class__(self.shape)
            for k, v in self.__data.items():
                out[k] = v
            for k, v in other.__data.items():
                out[k] = out[k] + v
        else: 
            raise Exception("Dimensions do not match.")
        

    
        