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

from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from itertools import product, combinations


fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_aspect("equal")

# draw cube
r = [-1, 1]
for s, e in combinations(np.array(list(product(r, r, r))), 2):
    if np.sum(np.abs(s-e)) == r[1]-r[0]:
        ax.plot3D(*zip(s, e), color="b")

# draw sphere
u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
x = np.cos(u)*np.sin(v)
y = np.sin(u)*np.sin(v)
z = np.cos(v)
ax.plot_wireframe(x, y, z, color="r")

# draw a point
ax.scatter([0], [0], [0], color="g", s=100)

# draw a vector
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d


class Arrow3D(FancyArrowPatch):

    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        FancyArrowPatch.draw(self, renderer)

a = Arrow3D([0, 1], [0, 1], [0, 1], mutation_scale=20,
            lw=1, arrowstyle="-|>", color="k")
ax.add_artist(a)
plt.show()