import numpy as np


EPS = 1e-15

class Vec3(object):
    # static methods manipulate a array
    @staticmethod
    def are_equal(a, b):
        return all(np.abs(a - b) <= EPS)

    @staticmethod
    def are_parallel(a, b):
        return all(np.cross(a, b) <= EPS)
    @staticmethod
    def are_perpendicular(a, b):
        return a.dot(b) <= EPS

    @staticmethod
    def normalize(v):
        # print("vec3::normalize::v:", v)
        if not np.any(v):
            return np.array(v)
        return np.array(v) / np.linalg.norm(v)

    @staticmethod
    def distance(p1, p2):
        np.sqrt(p2 - p1)
