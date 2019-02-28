"""
case one: no interfaces, ray casted from one transducer to infinite markoil
medium sample on cube and compare result from traditional method
"""

import getopt
import sys
import time
from multiprocessing import (Lock, Manager, Pool, Process, Queue, Value,
                             current_process, get_logger, log_to_stderr)

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from pylint.lint import multiprocessing

from pyHIFU import HIFU
from pyHIFU.box import Box, Sparse3D
from pyHIFU.io.config import readjson
from pyHIFU.physics.medium import InitMedium, MediaComplex
from pyHIFU.transducer import TElement, Transducer
from pyHIFU.visualization.figure import create_ax
from pyHIFU.visualization.mkplots import (plot_boundary, plot_box,
                                          plot_transducer)
from pyHIFU.visualization.plot_tensor import (plot_pressure_surf,
                                              plot_sliced_tensor)


def run_hifu(json_path, pyd_path, verbose=False, n_core=4):
    start_time = time.time()
    config = readjson(json_path=json_path)

    transducer_config = config['transducer']
    coor_file_path = "data/transducer_position.txt"
    textarray = np.loadtxt(coor_file_path)
    coordinates = textarray[::, 1:]
    T = Transducer(element_coordinates=coordinates, **transducer_config)

    mc = MediaComplex(config_json=config)

    init_medium_config = config['init_medium']
    init_medium = InitMedium.new_markoil(init_medium_config['boundary'])

    myhifu = HIFU(T)

    box_config = config['box']
    B = Box(*box_config['min'], *box_config['max'], box_config['step'])
    esp = transducer_config["element_init_paras"]
    cp = myhifu.run(
        init_medium,
        mc,
        B,
        esp['n_rays'],
        esp['trident_angle'],
        esp['theta_max'],
        n_core=n_core)
    plot_sliced_tensor(np.abs(cp), slicing_axis=2)


if __name__ == "__main__":
    config_path = 'data/case1.json'
    pyd_path = None

    options, remainder = getopt.getopt(sys.argv[1:], 'i:o:',
                                       ['output=', 'input='])
    for opt, arg in options:
        if opt in ('-o', '--output'):
            pyd_path = arg
        elif opt in ('-i', '--input'):
            config_path = arg

    print('using', config_path)
    run_hifu(config_path, pyd_path, verbose=False)
