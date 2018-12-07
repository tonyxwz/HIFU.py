import numpy as np


EPS = np.finfo(float).eps

class Vec3(object):
    # static methods manipulate a array
    @staticmethod
    def are_equal(a, b):
        return all(np.abs(a - b) < np.finfo(float).eps)

    @staticmethod
    def are_parallel(a, b):
        return all(np.cross(a, b) == 0)
    
    @staticmethod
    def are_perpendicular(a, b):
        return a.dot(b) == 0

    @staticmethod
    def normalize(v):
        return np.array(v) / np.linalg.norm(v)
        
    @staticmethod
    def distance(p1, p2):
        np.sqrt(p2 - p1)