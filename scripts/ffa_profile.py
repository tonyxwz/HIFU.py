import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pyHIFU.transducer import Transducer, TElement
from pyHIFU.io.config import readjson
from pyHIFU.physics.medium import Medium

config = readjson(json_path='../data/case1.json')

transducer_config = config['transducer']
coor_file_path = "../data/transducer_position.txt"
textarray = np.loadtxt(coor_file_path)
coordinates = textarray[::, 1:]

init_medium_config = config['medium_list'][0]
init_medium = Medium(**init_medium_config)

te = TElement(0, [0,0,0], 3.5e-3, 40, 1.2e6, np.array([0.14, 0,0]))
te.initialize(init_medium, 0, n=1, trident_angle=5e-2, theta_max=np.pi/20)

delta = 0.001
x = np.arange(0.02, 0.15, delta)
y = np.arange(-0.03, 0.03, delta)

X, Y = np.meshgrid(x, y)

z = np.zeros(X.shape, dtype=np.complex128)
for i in range(len(x)):
    for j in range(len(y)):
        coor = [x[i], y[j], 0]
        pc = te.ffa(coor)
        z[j, i] = pc

Z = np.abs(z)

fig, ax = plt.subplots()
CS = ax.contour(X, Y, Z)
ax.clabel(CS, inline=1, fontsize=10)
ax.set_title('absolute pressure contour')
plt.show()
