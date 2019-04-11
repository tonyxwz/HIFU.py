from pyHIFU.ray import Trident
import numpy as np
import matplotlib.pyplot as plt
from pyHIFU.physics.medium import MediaComplex
from pyHIFU.geometric.vec3 import Vec3
from pyHIFU.geometric.lines import Ray
from pyHIFU.physics import LONGITUDINAL, SHEAR
from pyHIFU.box import Box
from pyHIFU.box.sparse3d import Sparse3D
from pyHIFU.visualization.mkplots import plot_trident
from pyHIFU.visualization.plot_tensor import plot_sliced_tensor


medium_list = [
    {
        "is_init": True,
        "material_name": "markoil",

        "geometry": {
            "shape_type": "cuboid",
            "parameters": {
                "o": [0, -1, -1],
                "a": [0.1, 0, 0],
                "b": [0, 2, 0],
                "c": [0, 0, 2]
            }
        }
    },
    {
        "is_init": False,
        "material_name": "fake_bone",

        "geometry": {
            "shape_type": "cuboid",
            "parameters": {
                "o": [0.1, -1, -1],
                "a": [0.06, 0, 0],
                "b": [0, 2, 0],
                "c": [0, 0, 2]
            }
        }
    }
]

if __name__ == "__main__":
    # (x, y, z)
    start_coord = np.array([0, 0.01, 0])
    focus_coord = np.array([0.14, 0, 0])
    trident_angle = 5e-3
    pow_dire = Vec3.normalize(focus_coord - start_coord)

    ray_helper = Ray(start_coord, pow_dire)
    n_helper = ray_helper.perpendicularDirection()
    n2_helper = np.cross(pow_dire, n_helper)
    l = np.tan(trident_angle) * np.linalg.norm(pow_dire)
    assert l > 0
    p_end = start_coord + pow_dire
    a1_end = p_end + Vec3.normalize(n_helper) * l
    a2_end = p_end + Vec3.normalize(n2_helper) * l

    mc = MediaComplex(medium_list=medium_list)
    trident = Trident(
        start_coord,
        pow_dire,
        start_coord,
        a1_end - start_coord,
        start_coord,
        a2_end - start_coord,
        1,  # I0 = 1
        1.2e6,  # frequency
        0,  # initial phase
        len0=0.01,
        medium=mc[0],
        legacy=[],
        wave_type=LONGITUDINAL
    )
    tr_list_reflect = trident.reflect(mc)
    tr_list_refract = trident.refract(mc)

    assert tr_list_refract[0].pow_ray.d[2] == 0  # all intersection happen in x-y plane
    assert len(tr_list_reflect) == 0  # no ray reflected in markoil

    # print(trident.pow_ray.end, trident.bundle_identifier)
    # print(tr_list_refract[0].pow_ray.start, tr_list_refract[0].bundle_identifier)

    # box_min = trident.pow_ray.end - 0.01
    # box_max = trident.pow_ray.end + 0.01
    # l = .0001

    # box = Box(*box_min, *box_max, l=l)
    #
    # pc = np.zeros(box.nxyz, dtype=np.complex128)
    #
    # I = Sparse3D(box.nxyz)
    # ph = Sparse3D(box.nxyz)
    # counter = Sparse3D(box.nxyz, dtype=int)
    # box.intersect_trident(trident, I, ph, counter)
    # for k in I.getdata():
    #     # one could just assume I, ph, counter always have the same keys
    #     # use the Z of last tr because one bundle have the same medium
    #     pc[k] += np.sqrt(2 * trident.medium.Z * I[k] / counter[k]) * np.exp(
    #         1j * ph[k] / counter[k])
    #
    # I = Sparse3D(box.nxyz)
    # ph = Sparse3D(box.nxyz)
    # counter = Sparse3D(box.nxyz, dtype=int)
    # box.intersect_trident(tr_list_refract[0], I, ph, counter)
    # for k in I.getdata():
    #     # one could just assume I, ph, counter always have the same keys
    #     # use the Z of last tr because one bundle have the same medium
    #     pc[k] += np.sqrt(2 * tr_list_refract[0].medium.Z * I[k] / counter[k]) * np.exp(
    #         1j * ph[k] / counter[k])
    #
    # pressure = np.abs(pc)
    # print(pc.max())

    intensity_list = list()
    xs1 = np.arange(trident.pow_ray.endt - 0.01, trident.pow_ray.endt, 0.001)
    for x in xs1:
        intensity_list.append(trident.get_intensity_at(x))

    xs2 = np.arange(0, 0.01, 0.001)
    for x in xs2:
        intensity_list.append(tr_list_refract[0].get_intensity_at(x))

    xs = np.concatenate((xs1, xs2+trident.pow_ray.endt))

    fig = plt.figure()
    ax = fig.add_subplot(121, projection='3d')
    ax.set_aspect("equal")
    ax.set_xlim(0, 0.14)
    ax.set_ylim(-0.07, 0.07)
    ax.set_zlim(-0.07, 0.07)
    plot_trident(trident, ax)
    plot_trident(tr_list_refract[0], ax)

    ax2 = fig.add_subplot(122)
    ax2.plot(xs, intensity_list)
    plt.show()

