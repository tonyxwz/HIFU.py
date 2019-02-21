import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm


def plot_pressure_surf(t, i, ax):

    x = np.arange(t.shape[0])
    y = np.arange(t.shape[1])

    X, Y = np.meshgrid(y, x)
    Z = t[:,:,i]

    surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,linewidth=0, antialiased=False)
    # fig.colorbar(surf, shrink=0.5, aspect=5)
    return surf

if __name__ == "__main__":
    t = np.load("pressure_l0.0002_n500_t1550743618.npy")
    from pyHIFU.visualization.figure import create_ax
    fig = plt.figure()
    ax = create_ax(fig)
    surf = plot_pressure_surf(t, 15, ax)
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()
