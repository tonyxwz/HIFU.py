import numpy as np


EPS = 1e-10

class Vec3(object):
    # static methods manipulate a array
    @staticmethod
    def are_equal(a, b):
        assert len(a) == len(b) == 3
        return all(np.abs(a - b) <= EPS)
    @staticmethod
    def are_parallel(a, b):
        assert len(a) == len(b) == 3
        return all(np.cross(a, b) <= EPS)
    @staticmethod
    def are_perpendicular(a, b):
        assert len(a) == len(b) == 3
        return a.dot(b) <= EPS

    @staticmethod
    def normalize(v):
        # print("vec3::normalize::v:", v)
        assert len(v) == 3
        if not np.any(v):
            return np.array(v)
        return np.array(v) / np.linalg.norm(v)

    @staticmethod
    def distance(p1, p2):
        assert len(p1) == len(p2) == 3
        np.sqrt(p2 - p1)
    @staticmethod
    def is_infinitesimal(n):
        return n < EPS

    @staticmethod
    def is_zero(n):
        return Vec3.is_infinitesimal(n)
    @staticmethod
    def greaterthan(a, b):
        return a - b > EPS
    @staticmethod
    def num_are_equal(a, b):
        return np.abs(a - b) <= EPS
    @staticmethod
    def rotation_matrix(axis, theta):
        """
        Return the rotation matrix associated with counterclockwise rotation about
        the given axis by theta radians.
        """
        axis = np.asarray(axis)
        axis = axis / np.sqrt(np.dot(axis, axis))
        a = np.cos(theta / 2.0)
        b, c, d = -axis * np.sin(theta / 2.0)
        aa, bb, cc, dd = a * a, b * b, c * c, d * d
        bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
        return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                        [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                        [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

    @staticmethod
    def rotate(v, axis, theta):
        return np.dot(Vec3.rotation_matrix(axis, theta), v)

    @staticmethod
    def alt_rotation_matrix(axis, theta):
        from scipy.linalg import expm, norm
        from numpy import cross, eye, dot
        return expm(cross(eye(3), axis/norm(axis)*theta))

    @staticmethod
    def angle_between(v1, v2):
        assert len(v1) == len(v2) == 3
        v1 = Vec3.normalize(v1)
        v2 = Vec3.normalize(v2)
        costheta = np.dot(v1, v2)
        return np.arccos(costheta)
