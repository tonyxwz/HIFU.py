import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from matplotlib import cm
import numpy as np


def plot_sliced_tensor(t, slicing_axis=0):
    """
    `t`: 3D tensor data
    """
    axis = ['x', 'y', 'z']
    n = t.shape[slicing_axis]
    step = 1
    init = int(n / 2)
    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.2)
    im = plt.imshow(t[init,:,:])
    if slicing_axis == 1:
        im = plt.imshow(t[:,init,:])
    elif slicing_axis == 2:
        im = plt.imshow(t[:,:,init])
    plt.axis('off')

    ax_s = plt.axes([0.15, 0.07, 0.7, 0.03], facecolor="lightgoldenrodyellow")
    slider = Slider(ax_s, axis[slicing_axis]+'-Slice', 1, n, valinit=init, valstep=step, valfmt="%d")
    update_func = lambda val: update(val, t, im, slicing_axis=slicing_axis)
    slider.on_changed(update_func)
    plt.show()

def update(val, t, im, slicing_axis=0):
    new_slice = t[int(val-1),:,:]
    if slicing_axis == 1:
        new_slice = t[:,int(val-1),:]
    elif slicing_axis == 2:
        new_slice = t[:,:,int(val-1)]
    im.set_data(new_slice)


def plot_pressure_surf(t, ax, i=None):

    x = np.arange(t.shape[0])
    y = np.arange(t.shape[1])

    X, Y = np.meshgrid(y, x)
    if i is None:
        i = int(t.shape[2] / 2)
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


if __name__ == "__main__":
    from skimage import io
    vol = io.imread("https://s3.amazonaws.com/assets.datacamp.com/blog_assets/attention-mri.tif")
    plot_sliced_tensor(vol, slicing_axis=2)