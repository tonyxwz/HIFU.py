import numpy as np


class Point(np.array):
    def __init__(self, p):
        super().__init__(p)
        