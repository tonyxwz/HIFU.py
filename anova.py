"""
ANOVA one: no interfaces, ray casted from one transducer to infinite markoil
medium sample on cube and compare result from traditional method
"""

import getopt
import sys
import time

import numpy as np

from pyHIFU import HIFU
from pyHIFU.box import Box, Sparse3D
from pyHIFU.io.config import readjson, writejson
from pyHIFU.physics.medium import MediaComplex, Medium
from pyHIFU.transducer import TElement, Transducer

from pyHIFU.visualization.plot_tensor import (plot_pressure_surf,
                                              plot_sliced_tensor)
from evaluation import pnorm_distance
# try:
#     import seaborn
#     seaborn.set()
# except:
#     pass


def run_hifu(config_path, pyd_path, verbose=False, n_core=4, n_rays=None, trident_angle=None, theta_max=None):
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
    n_rays = n_rays or esp['n_rays']
    trident_angle = trident_angle or esp['trident_angle']
    theta_max = theta_max or esp['theta_max']
    cp = myhifu.run(
        mc,
        B,
        n_rays,
        trident_angle,
        theta_max,
        n_core=n_core)
    end_time = time.time()

    sig = "_".join([str(n_rays), str(trident_angle), str(theta_max)])
    print(f'use {end_time - start_time} seconds')
    np.save(f'npydata/ANOVA/complex_pressure_{start_time}_{sig}', cp)
    # plot_sliced_tensor(np.abs(cp), slicing_axis=2)
    return cp


if __name__ == "__main__":
    # import glob
    config_path = 'data/case1.json'
    pyd_path = None
    n_core = 8
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
    ground_truth = np.load(r"C:\Users\20175506\OneDrive - TU Eindhoven\Documents\TUE\1nternship\HIFU-ModenaEtAl-python\pressure_daniela_20000_markoil.npy")
    ground_truth = np.abs(ground_truth)
    ground_truth = ground_truth[1:, 1:, 1:]

    d_v_n_rays = dict()
    for n_rays in range(500, 20000, 1000):
        try:
            # npypath = glob.glob('npydata/ANOVA/complex_pressure_*_'+str(n_rays)+"_*.npy")
            # complexpressure = np.load(npypath[0])
            complexpressure = run_hifu(config_path, pyd_path, verbose=False, n_core=n_core, n_rays=n_rays)
            tmp_pressure = np.abs(complexpressure)
            tmp_dist = pnorm_distance(tmp_pressure, ground_truth)
            d_v_n_rays[n_rays] = tmp_dist
        except Exception as e:
            print(f"ERROR, n_rays={n_rays}, {e}")
    writejson(d_v_n_rays, "data/d_v_n_rays.json")

    d_v_trident_angle = dict()
    for trident_angle in np.arange(1e-3, 1e-2, 1e-3):
        try:
            complexpressure = run_hifu(config_path, pyd_path, verbose=False, n_core=n_core, trident_angle=trident_angle)
            tmp_pressure = np.abs(complexpressure)
            tmp_dist = pnorm_distance(tmp_pressure, ground_truth)
            d_v_trident_angle[trident_angle] = tmp_dist
        except Exception as e:
            print(f"ERROR, trident_angle={trident_angle}, {e}")
    writejson(d_v_trident_angle, "data/d_v_trident_angle.json")

    d_v_theta_max = dict()
    for theta_max in np.arange(1e-2, 1e-1, 1e-2):
        try:
            complexpressure = run_hifu(config_path, pyd_path, verbose=False, n_core=n_core, theta_max=theta_max)
            tmp_pressure = np.abs(complexpressure)
            tmp_dist = pnorm_distance(tmp_pressure, ground_truth)
            d_v_theta_max[theta_max] = tmp_dist
        except Exception as e:
            print(f"ERROR, theta_max={theta_max}, {e}")
    writejson(d_v_theta_max, "data/d_v_theta_max.json")
