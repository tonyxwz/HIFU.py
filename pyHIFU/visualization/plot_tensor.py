import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.animation as animation
import numpy as np


def plot_sliced_tensor(t, slicing_axis=0, update_cbar=False, ax=None):
    """
    `t`: 3D tensor data
    """
    axis = ['x', 'y', 'z']
    n = t.shape[slicing_axis]
    step = 1
    init = int(n / 2)
    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.2)
    if slicing_axis == 0:
        im = plt.imshow(t[init,:,:])
    if slicing_axis == 1:
        im = plt.imshow(t[:,init,:])
    elif slicing_axis == 2:
        im = plt.imshow(t[:,:,init])
    cb = plt.colorbar(im)
    # TODO add axis to real scale, slicing also show real scale
    plt.axis('off')

    ax_s = plt.axes([0.15, 0.07, 0.7, 0.03], facecolor="lightgoldenrodyellow")
    slider = Slider(ax_s, axis[slicing_axis]+'-Slice', 1, n, valinit=init, valstep=step, valfmt="%d")
    update_func = lambda val: update_im(val, t, im, cb, slicing_axis=slicing_axis, update_cbar=update_cbar)
    slider.on_changed(update_func)
    plt.show()


def update_im(val, t, im, cb, slicing_axis=0, update_cbar=False):
    if slicing_axis == 0:
        new_slice = t[int(val-1),:,:]
    if slicing_axis == 1:
        new_slice = t[:,int(val-1),:]
    elif slicing_axis == 2:
        new_slice = t[:,:,int(val-1)]
    im.set_data(new_slice)
    if update_cbar:
        cmax = new_slice.max()
        cmin = new_slice.min()
        cb.set_clim(cmin, cmax)
        cb_ticks = np.linspace(cmin, cmax, num=20, endpoint=True)
        cb.set_ticks(cb_ticks)
        cb.draw_all()


def plot_pressure_surf(Z):

    x = np.arange(Z.shape[0])
    y = np.arange(Z.shape[1])
    iz = int(Z.shape[2]/2)
    X, Y = np.meshgrid(y, x)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(X, Y, Z[:,:,iz], cmap=cm.coolwarm,linewidth=0, antialiased=False)
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()




if __name__ == "__main__":
    cp = np.load(r'D:/HIFU.py/npydata/complex_pressure_1551631771.9010646.npy')
    t = np.abs(cp)
    plot_pressure_surf(t)
