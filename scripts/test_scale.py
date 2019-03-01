"""
check the correct scale for the resulted pressure data
"""

import numpy as np
from pyHIFU.visualization.plot_tensor import plot_sliced_tensor


# runpy('case1')
t = np.load('npydata/pressure_l0.0002_n2000_t1550839895.npy')
tmax =  (t**2).max()

