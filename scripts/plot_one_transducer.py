import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pyHIFU.transducer import Transducer, TElement
from pyHIFU.io.config import readjson
from pyHIFU.physics.medium import InitMedium

config = readjson(json_path='HIFU.py/data/test_case1.json')

transducer_config = config['transducer']
coor_file_path = "HIFU.py/data/transducer_position.txt"
textarray = np.loadtxt(coor_file_path)
coordinates = textarray[::,1:]

init_medium_config = config['init_medium']
init_medium = InitMedium.new_markoil(init_medium_config['boundary'])

te = TElement(0, [0,0,0], 3.5e-3, 40, 1.2e6, np.array([0.14, 0,0]))
te.initialize(init_medium, 0, n=3000, trident_angle=5e-2, theta_max=np.pi/20)

from pyHIFU.visualization.mkplots import plot_TElements
from pyHIFU.visualization.figure import create_ax

fig = plt.figure()
ax = create_ax(fig)

plot_TElements(te, ax)


plt.show()
