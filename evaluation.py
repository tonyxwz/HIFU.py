"""
Compare the similarity between the simulated result and Daniela's result
"""
import numpy as np


def pnorm_distance(a1, a2):
    assert a1.shape == a2.shape
    a1 = normalize(a1)
    a2 = normalize(a2)
    x = np.abs(a1 - a2)
    x = np.power(x, 3)
    s = np.sum(x)
    s = np.power(s, 1/3)
    return s


def normalize(a):
    m = a.min()
    M = a.max()
    a2 = (a - m) / (M - m)
    return a2
