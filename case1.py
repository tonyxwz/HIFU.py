# case one: no interfaces, ray casted from one transducer to cube medium
# compare result from traditional method

from pyHIFU.transducer import Transducer
from pyHIFU.io.config import readjson

config = readjson(json_path='data/example2.json')
transducer_config = config['transducer']
T = Transducer(transducer_config)

from pyHIFU.physics.medium import MediaComplex

mc = MediaComplex(config_json=config)

from math import inf
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
t = inf
interface = None

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')

n = 100
pps = []
pa1s = []
pa2s = []

for i,el in enumerate(T):
    print("element", i)
    for j,tr in enumerate(el):
        print("trident", j)
        for k, m in enumerate(mc):
            print("medium", k)
            for x,f in enumerate(m.shape):
                print("face:", x)
                pp = tr.pow_ray.intersect_plane(f)
                pa1 = tr.aux_ray1.intersect_plane(f)
                pa2 = tr.aux_ray2.intersect_plane(f)
                if (pp is not None):
                    l = np.linalg.norm(pp-tr.pow_ray.start)
                    if t > l:
                        t = l
                        interface = f
                    pps.append(pp)
                    print("pp:", pp)
                if pa1 is not None:
                    pa1s.append(pa1)
                    print("pa1", pa1)
                if pa2 is not None: 
                    pa2s.append(pa2)
                    print("pa2", pa2)
                    