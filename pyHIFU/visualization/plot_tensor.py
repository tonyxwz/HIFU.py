import matplotlib.pyplot as plt
from matplotlib.widgets import Slider


def plot_sliced_tensor(t, slicing_axis=0):
    """
    `t`: 3D tensor data
    """
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
    slider = Slider(ax_s, 'x-Slice', 1, n, valinit=init, valstep=step, valfmt="%d")
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


if __name__ == "__main__":
    from skimage import io
    vol = io.imread("https://s3.amazonaws.com/assets.datacamp.com/blog_assets/attention-mri.tif")
    plot_sliced_tensor(vol)