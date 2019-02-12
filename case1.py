"""
case one: no interfaces, ray casted from one transducer to infinite markoil
medium sample on cube and compare result from traditional method
"""
#%% 
import time
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from pyHIFU.transducer import Transducer
from pyHIFU.io.config import readjson
from pyHIFU.physics.medium import MediaComplex
from pyHIFU.physics.medium import InitMedium
from pyHIFU.box import Box
from pyHIFU.box import Sparse3D

from pyHIFU.visualization.mkplots import plot_transducer, plot_boundary, plot_box
from pyHIFU.visualization.figure import create_ax
from pyHIFU.visualization.plot_tensor import plot_sliced_tensor


#%% 
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

bundle_dict = T.cast()

end_time = time.time()
print("initialization and casting time:", end_time-start_time, "seconds")

#%%
# start sampling
start_time = time.time()
box_config = config['box']

B = Box(*box_config['min'], *box_config['max'], box_config['step'])
print(len(B.lattrix))
end_time = time.time()


def measure(box:Box, bd):
    pc = np.zeros(box.nxyz, dtype=np.complex128)  # complex pressure

    for _bundle_str, tr_list in bd.items():
        I = Sparse3D(box.nxyz)
        ph = Sparse3D(box.nxyz)
        counter = Sparse3D(box.nxyz, dtype=int)
        for tr in tr_list:
            box.intersect_trident(tr, I, ph, counter)

        # TODO update pc
        for k in I.getdata():
            # one could just assume I, ph, counter always have the same keys
            # use the Z of last tr because one bundle have the same medium
            pc[k] = np.sqrt(2 * tr.medium.Z * I[k] / counter[k]) * np.exp(1j * ph[k])

    return pc
            
complex_pressure = measure(B, bundle_dict)

print("sampling time:", end_time-start_time, "seconds")

real_pressure = np.abs(complex_pressure)
np.save('pressure', real_pressure)

plot_sliced_tensor(real_pressure)

verbose = False
if verbose:
    fig = plt.figure()
    ax = create_ax(fig, 121)
    plot_transducer(T, ax)
    plot_boundary(T.init_medium.boundary, ax)
    plot_box(B, ax, title="Box with Transducer")

    ax = create_ax(fig, 122)

    plt.show()
