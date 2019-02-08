# case one: no interfaces, ray casted from one transducer to infinite markoil medium
# sample on cube and compare result from traditional method

from pyHIFU.transducer import Transducer
from pyHIFU.io.config import readjson
import numpy as np
from pyHIFU.physics.medium import MediaComplex
from pyHIFU.physics.medium import InitMedium
from pyHIFU.box import Box

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pyHIFU.visualization.mkplots import plot_transducer, plot_boundary, plot_box
from pyHIFU.visualization.figure import create_ax
import time


start_time = time.time()
config = readjson(json_path='data/test_case1.json')

transducer_config = config['transducer']
coor_file_path = "data/transducer_position.txt"
textarray = np.loadtxt(coor_file_path)
coordinates = textarray[::,1:]
T = Transducer(element_coordinates=coordinates, **transducer_config)

mc = MediaComplex(config_json=config)

init_medium_config = config['init_medium']
init_medium = InitMedium.new_markoil(init_medium_config['boundary'])
T.initialize(init_medium, **transducer_config["element_init_paras"])

trident_dict = T.cast()

end_time = time.time()
print("initialization and casting time:", end_time-start_time, "seconds")


# start sampling
start_time = time.time()
B = Box(0.05,-0.025,-0.025, 0.1, 0.025, 0.025, 0.0005)
print(len(B.lattrix))
end_time = time.time()
print("sampling time:", end_time-start_time, "seconds")

fig = plt.figure()
ax = create_ax(fig)
plot_transducer(T, ax)
plot_boundary(T.init_medium.boundary, ax)
plot_box(B, ax, title="Box with Transducer")

plt.show()
