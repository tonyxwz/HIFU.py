# case one: no interfaces, ray casted from one transducer to infinite markoil medium
# sample on cube and compare result from traditional method

from pyHIFU.transducer import Transducer
from pyHIFU.io.config import readjson
import numpy as np
from pyHIFU.physics.medium import MediaComplex
from pyHIFU.box import Box

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pyHIFU.visualization.mkplots import *
import time


start_time = time.time()
config = readjson(json_path='data/test_case1.json')
transducer_config = config['transducer']
coor_file_path = "data/transducer_position.txt"
textarray = np.loadtxt(coor_file_path)
coordinates = textarray[::,1:]

mc = MediaComplex(config_json=config)
T = Transducer(mc[0], element_coordinates=coordinates, **transducer_config)
B = Box(0.1,-0.15,-0.15, 0.4, 0.15, 0.15, 0.02)
# B = Box(3, -2, -2, 7, 2, 2, 0.2)
T.initialize()

end_time = time.time()
print("Total time:", end_time-start_time, "seconds")


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
plot_transducer(T, ax)
plot_boundary(T.init_medium.boundary, ax)
plot_box(B, ax, title="Box with Transducer")

plt.show()
