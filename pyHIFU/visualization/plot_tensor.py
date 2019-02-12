import matplotlib.pyplot as plt
from matplotlib.widgets import Slider


def plot_sliced_tensor(t, ):
    """
    `t`: 3D tensor data
    """
    nx, ny, nz = t.shape
    step_z = 1
    init_z = int(nz / 2)
    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.2)
    im = plt.imshow(t[:,:,init_z])
    plt.axis('off')

    ax_s = plt.axes([0.15, 0.07, 0.7, 0.03], facecolor="lightgoldenrodyellow")
    slider = Slider(ax_s, 'z-Slice', 1, nz, valinit=init_z, valstep=step_z, valfmt="%d")
    update_func = lambda val: update(val, t, im)
    slider.on_changed(update_func)
    plt.show()

def update(val, t, im):
    new_slice = t[:,:,int(val-1)]
    im.set_data(new_slice)


if __name__ == "__main__":
    from skimage import io
    vol = io.imread("https://s3.amazonaws.com/assets.datacamp.com/blog_assets/attention-mri.tif")
    plot_sliced_tensor(vol)