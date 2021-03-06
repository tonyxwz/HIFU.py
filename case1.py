"""
case one: no interfaces, ray casted from one transducer to infinite markoil
medium sample on cube and compare result from traditional method
"""

import getopt
import sys
import time

import numpy as np

from pyHIFU import HIFU
from pyHIFU.box import Box, Sparse3D
from pyHIFU.io.config import readjson
from pyHIFU.physics.medium import MediaComplex, Medium
from pyHIFU.transducer import TElement, Transducer

from pyHIFU.visualization.plot_tensor import (plot_pressure_surf,
                                              plot_sliced_tensor)
# try:
#     import seaborn
#     seaborn.set()
# except:
#     pass


def run_hifu(config_path, pyd_path, verbose=False, n_core=4):
    start_time = time.time()
    config = readjson(json_path=config_path)

    transducer_config = config['transducer']
    coor_file_path = "data/transducer_position.txt"
    textarray = np.loadtxt(coor_file_path)
    coordinates = textarray[::, 1:]
    # coordinates = coordinates[9]
    # coordinates = np.reshape(coordinates, (1, 3))
    T = Transducer(element_coordinates=coordinates, **transducer_config)

    mc = MediaComplex(medium_list=config["medium_list"])

    box_config = config['box']
    B = Box(*box_config['min'], *box_config['max'], box_config['step'])

    esp = transducer_config["element_init_paras"]
    myhifu = HIFU(T)
    cp = myhifu.run(
        mc,
        B,
        esp['n_rays'],
        esp['trident_angle'],
        esp['theta_max'],
        n_core=n_core)
    end_time = time.time()
    print(f'use {end_time - start_time} seconds')
    np.save(f'npydata/complex_pressure_{start_time}', cp)
    plot_sliced_tensor(np.abs(cp), slicing_axis=2)


if __name__ == "__main__":
    config_path = 'data/case1.json'
    pyd_path = None
    n_core = 0
    options, remainder = getopt.getopt(sys.argv[1:], 'i:o:j:',
                                       ['output=', 'input=', 'n_core='])
    for opt, arg in options:
        if opt in ('-o', '--output'):
            pyd_path = arg
        elif opt in ('-i', '--input'):
            config_path = arg
        elif opt in ('-j', '--n_core'):
            n_core = int(arg)
    print('using', config_path)
    run_hifu(config_path, pyd_path, verbose=False, n_core=n_core)
