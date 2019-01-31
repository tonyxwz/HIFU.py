import numpy as np
from pyHIFU.geometric.surfaces import Plane
from pyHIFU.physics.medium import Medium
from pyHIFU.geometric.vec3 import Vec3
from pyHIFU.geometric.lines import Ray as GeoRay
# from numpy import linalg as LA


def snell(v_in, n, c1, c2, return_everything=False):
    """ Snell's law
    `ray`: pyHIFU.ray.PowRay / AuxRay
    `interface`: pyHIFU.geometric.surfaces.Plane
    `c1`, `c2`: speed
    return: list of geo rays, leave the rest of the job to the calling function
    """
    # # check the state of the two interfaces
    
    # # if s1 is solid -> 2 reflected rays (long and shear)
    # ip = ray.intersect_plane(interface) # intersection point as the starting of new rays
    if np.dot(v_in, n) > 0:  # pointing into medium1
        n = -n
    nbr = c2 / c1
    a = v_in - (np.dot(v_in, n))*n
    sinalpha_in = np.linalg.norm(a)
    sinalpha_out = nbr * sinalpha_in
    if return_everything:
        return sinalpha_out, a/sinalpha_in
    else:
        return sinalpha_out


def vray_reflected(ray, interface: Plane, c1=None, c2=None):
    """ vector reflected
    ray: incident ray
    interface: plane
    c1: sound speed for incident ray
    c2: sound speed for reflected ray
    """
    n = interface.n
    if np.dot(ray.d, n)>0:
        n = -n
    if c1 is None and c2 is None:
        # ordinary reflection
        v_out = ray.d - 2 * (np.dot(ray.d, n)) * n
        return v_out
    else:
        sintheta_out, vh = snell(ray.d, n, c1, c2, return_everything=True)
        if sintheta_out > 1:
            return None  # reflection not possible
        else:
            cosalpha_out = np.sqrt(1 - sintheta_out ** 2)
            v_out = cosalpha_out * n + sintheta_out * vh
            return v_out

def vpol_reflected(ray, interface: Plane, c1=None, c2=None):
    pass

def fresnel(ray):
    pass