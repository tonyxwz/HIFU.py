import numpy as np
import matplotlib.pyplot as plt
import seaborn
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d as mp3d
from random import random

def plot_transducer(T, ax):
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    k = T.element_coordinates

    X = k[:,0]
    Y = k[:,1]
    Z = k[:,2]
    ax.scatter(X, Y, Z, s=1, color="r")
    ax.scatter(*T.nature_focus, color="teal")
    X = np.append(X, T.nature_focus[0])
    Y = np.append(Y, T.nature_focus[1])
    Z = np.append(Z, T.nature_focus[2])

    max_range = np.array([X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max() / 2.0

    mid_x = (X.max()+X.min()) * 0.5
    mid_y = (Y.max()+Y.min()) * 0.5
    mid_z = (Z.max()+Z.min()) * 0.5
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)
    
    for te in T:
        plot_TElements(te, ax)

    # plt.show()

def plot_TElements(te, ax):
    # fig = plt.figure()
    # ax = fig.gca(projection='3d')
    for tr in te:
        xyz = np.concatenate((tr.pow_ray.p, tr.pow_ray.end))
        xyz = xyz.reshape((2,3))
        ax.plot(xyz[:,0], xyz[:,1], xyz[:,2], linewidth=0.5, color="m")

        xyz = np.concatenate((tr.aux_ray1.p, tr.aux_ray1.end))
        xyz = xyz.reshape((2,3))
        ax.plot(xyz[:,0], xyz[:,1], xyz[:,2], '--', linewidth=0.5, color="g")
        
        xyz = np.concatenate((tr.aux_ray2.p, tr.aux_ray2.end))
        xyz = xyz.reshape((2,3))
        ax.plot(xyz[:,0], xyz[:,1], xyz[:,2], '--', linewidth=0.5, color="g")

def plot_boundary(pl, ax):
    x = pl.p[0]
    vertices = [(x, -0.1, -0.1),
        (x, -0.1, 0.1),
        (x, 0.1, 0.1),
        (x, 0.1, -0.1)
        ]
    alpha = 0.3
    face = mp3d.art3d.Poly3DCollection([vertices], alpha=alpha, linewidth=1)

    face.set_facecolor((0, 0, 1, alpha))
    ax.add_collection3d(face)

def plot_box(b, ax, title=""):
    alpha = 0.3
    
    for face in b:
        # face: rectangle
        vertices = list()
        for edge in face.edges:
            vertices.append(edge.start)
        face = mp3d.art3d.Poly3DCollection([vertices], alpha=alpha, linewidth=1)
        face.set_facecolor((random(), random(), random(), alpha))
        ax.add_collection3d(face)
    
    if len(title):
        ax.set_title(title)
        rangex = b.o2[0] - b.o1[0]
        rangey = b.o2[1] - b.o1[1]
        rangez = b.o2[2] - b.o1[2]

        ax.set_xlim(b.o1[0]-0.3*rangex, b.o1[0]+1.3*rangex)
        ax.set_xlabel("x (m)")
        ax.set_ylim(b.o1[1]-0.3*rangey, b.o1[1]+1.3*rangey)
        ax.set_ylabel("y (m)")
        ax.set_zlim(b.o1[2]-0.3*rangez, b.o1[2]+1.3*rangez)
        ax.set_zlabel("z (m)")


def plot_trident(tr, ax):
    xyz = np.concatenate((tr.pow_ray.p, tr.pow_ray.end))
    xyz = xyz.reshape((2,3))
    ax.plot(xyz[:,0], xyz[:,1], xyz[:,2], color="m")

    xyz = np.concatenate((tr.aux_ray1.p, tr.aux_ray1.end))
    xyz = xyz.reshape((2,3))
    ax.plot(xyz[:,0], xyz[:,1], xyz[:,2], '--', color="g")
    
    xyz = np.concatenate((tr.aux_ray2.p, tr.aux_ray2.end))
    xyz = xyz.reshape((2,3))
    ax.plot(xyz[:,0], xyz[:,1], xyz[:,2], '--', color="g")

def plot_reflection(tr, interface, ax):
    # TODO: the reflection progress of one transducer on one interface
    # also give resulted shear ray and longitudinal ray different color
    pass
