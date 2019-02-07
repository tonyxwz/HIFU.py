import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def create_ax(fig, coordinates=111):
    """ even unit length for all axis """
    ax = fig.add_subplot(coordinates, projection="3d")
    ax.set_aspect("equal")
    return ax