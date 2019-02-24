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

from pyHIFU.box import Box, Sparse3D
from pyHIFU.io.config import readjson
from pyHIFU.physics.medium import InitMedium, MediaComplex
from pyHIFU.transducer import TElement, Transducer
from pyHIFU.visualization.figure import create_ax
from pyHIFU.visualization.mkplots import (plot_boundary, plot_box,
                                          plot_transducer)
from pyHIFU.visualization.plot_tensor import (plot_pressure_surf,
                                              plot_sliced_tensor)


def run(json_path, pyd_path, verbose=False):
    start_time = time.time()
    config = readjson(json_path=json_path)

    transducer_config = config['transducer']
    coor_file_path = "data/transducer_position.txt"
    textarray = np.loadtxt(coor_file_path)
    coordinates = textarray[::,1:]
    T = Transducer(element_coordinates=coordinates, **transducer_config)

    mc = MediaComplex(config_json=config)

    init_medium_config = config['init_medium']
    init_medium = InitMedium.new_markoil(init_medium_config['boundary'])

    T.initialize(init_medium, **transducer_config["element_init_paras"], verbose=verbose)
    end_time = time.time()
    if verbose: print("--- initialization:", end_time-start_time, "seconds ---")
    
    start_time = end_time
    bundle_dict = T.cast()
    end_time = time.time()
    if verbose: print("--- casting time:", end_time-start_time, "seconds ---")
    
    # start sampling
    start_time = end_time
    box_config = config['box']

    B = Box(*box_config['min'], *box_config['max'], box_config['step'])
    # print(len(B.lattrix))

    complex_pressure = measure(B, bundle_dict, verbose=verbose)

    end_time = time.time()
    if verbose: print("--- sampling time:", end_time-start_time, "seconds ---")
    
    if pyd_path is None:
        pyd_path = 'pressure'+'_l'+str(B.lx)+"_n"+str(len(T[0])) + "_t"+str(int(time.time()))
    real_pressure = np.abs(complex_pressure)
    np.save('npydata/'+pyd_path, real_pressure)
    np.save('npydata/'+'c'+pyd_path, complex_pressure)

    if verbose:
        # TODO plot surface
        plot_sliced_tensor(real_pressure, slicing_axis=2)

def measure(box:Box, bd, verbose=False):
    pc = np.zeros(box.nxyz, dtype=np.complex128)  # complex pressure
    c = np.zeros(box.nxyz, dtype=int)

    for _bundle_str, tr_list in bd.items():
        if verbose: print("Bundle:", _bundle_str)
        I = Sparse3D(box.nxyz)
        ph = Sparse3D(box.nxyz)
        counter = Sparse3D(box.nxyz, dtype=int)
        for tr in tr_list:
            box.intersect_trident(tr, I, ph, counter)

        for k in I.getdata():
            # one could just assume I, ph, counter always have the same keys
            # use the Z of last tr because one bundle have the same medium
            pc[k] += np.sqrt(2 * tr.medium.Z * I[k] / counter[k]) * np.exp(1j * ph[k] / counter[k])

    return pc


def run_mp(json_path, pyd_path, verbose=False, n_core=10):
    start_time = time.time()
    config = readjson(json_path=json_path)

    transducer_config = config['transducer']
    coor_file_path = "data/transducer_position.txt"
    textarray = np.loadtxt(coor_file_path)
    coordinates = textarray[::,1:]
    T = Transducer(element_coordinates=coordinates, **transducer_config)

    mc = MediaComplex(config_json=config)

    init_medium_config = config['init_medium']
    init_medium = InitMedium.new_markoil(init_medium_config['boundary'])

    T.initialize(init_medium, **transducer_config["element_init_paras"], verbose=verbose)
    end_time = time.time()
    if verbose: print("--- initialization:", end_time-start_time, "seconds ---")
    
    # start sampling
    start_time = end_time
    box_config = config['box']
    B = Box(*box_config['min'], *box_config['max'], box_config['step'])

    # parallelization, give transducer casting and measuring job to process pool
    pool = Pool(n_core)  # process worker pool, 10 CPU cores
    m = Manager()
    cpq = m.Queue()  # "Q"ueue for process to communicate complex "P"ressure
    lock = Lock()

    for te in T:
        pool.apply_async(measure_kernel, args=(te, B, mc, cpq, verbose))
        # print(ar.get())
    pool.close()
    pool.join()
    if verbose: print("---  all process finished  ---")

    complex_pressure = np.zeros(B.nxyz, dtype=np.complex128)
    
    while not cpq.empty():
        complex_pressure += cpq.get()
    
    end_time = time.time()
    if verbose: print("--- casting time:", end_time-start_time, "seconds ---")
    
    # save ndarray to file
    if pyd_path is None:
        pyd_path = 'pressure'+'_l'+str(B.lx)+"_n"+str(len(T[0])) + "_t"+str(int(time.time()))
    real_pressure = np.abs(complex_pressure)
    np.save('npydata/'+pyd_path, real_pressure)
    np.save('npydata/'+'c'+pyd_path, complex_pressure)

    if verbose:
        # TODO plot surface
        plot_sliced_tensor(real_pressure, slicing_axis=2)


def measure_kernel(te: TElement, box: Box, mc: MediaComplex, q, verbose):
    """ measuring process kernel function """
    pc = np.zeros(box.nxyz, dtype=np.complex128)
    pname = current_process().name
    if verbose: print(f"TE #{te.el_id} is assigned to proc {pname}")
    # TODO move te.initialize to here
    bundle_dict = te.cast(mc)
    if verbose: print(f"TE #{te.el_id} casted")
    for _bundle_str, tr_list in bundle_dict.items():
        if verbose: print(f"Bundle{_bundle_str}")
        I = Sparse3D(box.nxyz)
        ph = Sparse3D(box.nxyz)
        counter = Sparse3D(box.nxyz, dtype=int)
        for tr in tr_list:
            box.intersect_trident(tr, I, ph, counter)
        
        for k in I.getdata():
            # one could just assume I, ph, counter always have the same keys
            # use the Z of last tr because one bundle have the same medium
            pc[k] += np.sqrt(2 * tr.medium.Z * I[k] / counter[k]) * np.exp(1j * ph[k] / counter[k])
    q.put(pc)
    

if __name__ == "__main__":
    config_path = 'data/test_case1.json'
    pyd_path = None

    options, remainder = getopt.getopt(sys.argv[1:], 'i:o:', ['output=', 'input='])
    for opt, arg in options:
        if opt in ('-o', '--output'):
            pyd_path = arg
        elif opt in ('-i', '--input'):
            config_path = arg
    print("500 rays per transducer for comparison")
    t1 = time.time()
    run(config_path, pyd_path, verbose=False)
    t2 = time.time()
    print(f"serial processing cost {t2-t1} seconds")

    t3 = time.time()
    run_mp(config_path, pyd_path, verbose=False, n_core=12)
    t4 = time.time()
    print(f"multi-processing cost {t4-t3} seconds")