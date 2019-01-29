# case one: no interfaces, ray casted from one transducer to infinite markoil medium
# sample on cube and compare result from traditional method

from pyHIFU.transducer import Transducer
from pyHIFU.io.config import readjson
import numpy as np
from pyHIFU.physics.medium import MediaComplex

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pyHIFU.visualization.mkplots import *


config = readjson(json_path='data/test_case1.json')
transducer_config = config['transducer']
mc = MediaComplex(config_json=config)
coor_file_path = "data/transducer_position.txt"
textarray = np.loadtxt(coor_file_path)
coordinates = textarray[::3,1:]

T = Transducer(mc[0], element_coordinates=coordinates, **transducer_config)
T.initialize()
T.initial_cast()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
plot_transducer(T, ax)
plot_boundary(T.init_medium.boundary, ax)

plt.show()




                    