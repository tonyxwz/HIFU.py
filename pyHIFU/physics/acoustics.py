from pyHIFU.geometric.surfaces import Plane
from pyHIFU.physics.medium import Medium
from pyHIFU.geometric.vec3 import Vec3
from pyHIFU.geometric.lines import Ray as GeoRay
# from numpy import linalg as LA
import numpy as np
from numpy import linalg as LA
import cmath as cmath


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
    a = v_in - (np.dot(v_in, n))*n  # horizontal vector length is sin(alpha_in)
    sinalpha_in = np.linalg.norm(a)
    sinalpha_out = nbr * sinalpha_in
    if return_everything:
        return sinalpha_out, a/sinalpha_in  # horizontal vector length is 1
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
    if np.dot(ray.d, n) > 0:
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
    """polarization direction of reflected ray (shear only)
    
    Parameters
    ----------
    ray : Ray
    interface : Plane
    c1 : sound speed of first medium, optional
    c2 : sound speed of second medium, optional    
    """

    pass


def vray_refracted(ray, interface: Plane, c1=None, c2=None):
    """ vector refracted
    ray: incident ray
    interface: plane
    c1: sound speed for incident ray
    c2: sound speed for reflected ray
    """
    n = interface.n
    if np.dot(ray.d, n) > 0:
        n = -n

    sintheta_out, vh = snell(ray.d, n, c1, c2, return_everything=True)
    if sintheta_out > 1:
        return None  # total internal reflection, no refraction
    else:
        cosalpha_out = np.sqrt(1 - sintheta_out ** 2)
        v_out = - cosalpha_out * n + sintheta_out * vh
        return v_out


def vpol_refracted(ray, interface: Plane, c1=None, c2=None):
    """polarization direction of refracted ray (shear only)

    Parameters
    ----------
    ray : Ray
    interface : Plane
    c1 : sound speed of first medium, optional
    c2 : sound speed of second medium, optional
    """
    pass


# -----------------------------------------
# | begin of Daniela's acoustic functions |
# -----------------------------------------
def coefficient_ll(v_in, n, c1, c2, rho1, rho2):
    """liquid to liquid
    """
    v_in = v_in / LA.norm(v_in)
    n = n / LA.norm(n)

    if np.dot(v_in, n) < 0:
        n = -n

    costheta_i = np.dot(v_in, n)
    sintheta_i = np.sqrt(1 - costheta_i**2)
    tmp = (c2 / c1) * sintheta_i
    if tmp > 1:
        # Total internal reflection
        sintheta_t = 1
        costheta_t = 0
    else:
        sintheta_t = tmp
        costheta_t = np.sqrt(1 - sintheta_t**2)

    T = (4 * rho1 * c1 * rho2 * c2 * costheta_t * costheta_i) / (rho2 * c2 * costheta_i + rho1 * c1 * costheta_t) ** 2
    R = ((rho2 * c2 * costheta_i - rho1 * c1 * costheta_t) ** 2) / (rho2 * c2 * costheta_i + rho1 * c1 * costheta_t) ** 2

    Tp = (4 * rho1 * c1 * rho2 * c2 * costheta_t * costheta_i) / (rho2 * c2 * costheta_i + rho1 * c1 * costheta_t)
    Rp = (rho2 * c2 * costheta_i - rho1 * c1 * costheta_t) / (rho2 * c2 * costheta_i + rho1 * c1 * costheta_t)
    # phase of reflection coeff:
    Ph_refl = np.angle(Rp)
    Ph_tr = np.angle(Tp)
    T = T * costheta_t / costheta_i
    return T, R, Ph_refl, Ph_tr


def refraction0(v_in, n, c_in, c_out):
    v_in = v_in / LA.norm(v_in)
    n = n / LA.norm(n)

    if np.dot(v_in, n) < 0:
        n = -n  #directed into out

    nbr = c_out / c_in
    a = v_in - (np.dot(v_in,n))*n
    na = LA.norm(a)
    if na == 0:  # incoming ray parallel to normal
        refr = True
        v_out = v_in
    else:  # incoming ray not perpendicular
        sintheta2 = nbr * na  # Snellius: sintheta1 / sintheta2 = c_in / c_out
        if sintheta2 <= 1:  # refraction
            costheta2 = np.sqrt(1 - sintheta2 ** 2)
            refr = True
            v_out = costheta2 * n + nbr * a
        else:  # No refraction
            refr = False
            v_out = np.array([0, 0, 0])

    return refr, v_out


def reflection(v_in, n):
    v_in = v_in / LA.norm(v_in)
    n = n / LA.norm(n)
    v_out = v_in - 2 * (np.dot(v_in,n))*n
    return v_out


def M2BCoef(alpha_in, rho_liquid, rho_solid, c_long_liquid, c_long_solid, c_shear_solid):
    # COMPUTATION OF REFLECTION AND TRANSMISSION COEF VAN LIQUID(LONGIT INPUT) TO
    # SOLID(LONG EN SHEAR OUTPUT) GIVES POWER COEFF: Refl, Transm_long, Transm_shear
    # AND VELOCITY (= DISPLACEMENT) COEFF VRefl, VTransm_long, VTransm_shear
    # alpha_in in radials!!
    # Lamé params solid:
    mu_solid = rho_solid * c_shear_solid ** 2
    lambda_solid = rho_solid * c_long_solid ** 2 - 2 * mu_solid
    # Lamé params liquid:
    lambda_liquid = rho_liquid * c_long_liquid ** 2
    # compute critical angles:
    alpha_long_crit = np.arcsin(c_long_liquid / c_long_solid)
    alpha_shear_crit = np.arcsin(c_long_liquid / c_shear_solid)
    alpha_long_crit_degr = alpha_long_crit * 180 / np.pi
    alpha_shear_crit_degr = alpha_shear_crit * 180 / np.pi

    # compute beta_long: Snellius: sin(alpha_in) / sin(beta_long) = c_long_liquid / c_long_solid;

    sinbeta_long = np.sin(alpha_in) * c_long_solid / c_long_liquid
    cosbeta_long = cmath.sqrt(1 - sinbeta_long ** 2) # pure imag if alpha_in > alpha_long_crit
    # compute beta_shear: Snellius: sin(alpha_in) / sin(beta_shear) = c_long_liquid / c_shear_solid;
    sinbeta_shear = np.sin(alpha_in) * c_shear_solid / c_long_liquid
    cosbeta_shear = cmath.sqrt(1 - sinbeta_shear ** 2) # pure imag if alpha_in > alpha_shear_crit

    # FORMULE VERSION(FROM pdf with unknown author)
    Z1 = rho_liquid * c_long_liquid / np.cos(alpha_in)
    ZL = rho_solid * c_long_solid / cosbeta_long
    Zs = rho_solid * c_shear_solid / cosbeta_shear
    Zeff = ZL * (cosbeta_shear ** 2 - sinbeta_shear ** 2) ** 2 + Zs * (2 * sinbeta_shear * cosbeta_shear) ** 2;

    # REFLECTION:
    Rdf = (Zeff - Z1) / (Zeff + Z1)
    # this is the reflection coeff for the VELOCITY POTENTIAL(see pdf) BUT THAT IS EQUAL TO THE
    # COEFF FOR DISPLACEMENTS OR VELOCITIES
    VRefl = Rdf
    # now compute power refl coeff(note: same inpedances and same in and  outgoing angles
    # at reflection!):
    Refl = np.abs(Rdf) ** 2

    # TRANSMISSION, LONGITUDINAL:
    Tdf = (rho_liquid / rho_solid) * 2 * ZL * (cosbeta_shear ** 2 - sinbeta_shear ** 2) / (Zeff + Z1)
    # Take care: this formula gives again the transm coeff for the VELOCITY  POTENTIAL(see
    # pdf).The matrix computation in Auld or Rose gives immediately the coeff for displacement and velocity.
    # BUT: Tdf, velocity(Auld) = Tdf, vel potential * | | k_shear_solid | | / | | k_long_liquid | |
    # = Tdf, vel potential * c_long_liquid / c_long_solid
    Tdf = Tdf * c_long_liquid / c_long_solid
    # now Tdf is the displacement( or velocity) transmission coeff, like in Auld, for longit waves
    VTransm_long = Tdf
    # now compute power transmission coeff for longit(using impedances and angles):
    Tdpow = np.abs(Tdf) ** 2 * (rho_solid * c_long_solid) / (rho_liquid * c_long_liquid) * cosbeta_long / np.cos(alpha_in)

    # Note: Tdpow is real below the critical angle.Above the critical angle is  cosbeta_long
    # pure imaginairy and hence also Tdpow is pure imaginairy. In that case there is no
    # power transmission(evanescent waves), so we can take the reall part of Tdpow to
    # describe the actual power transmission:
    Transm_long = np.real(Tdpow)

    # TRANSMISSION, SHEAR:
    Tsf = -(rho_liquid / rho_solid) * 2 * Zs * (2 * sinbeta_shear * cosbeta_shear) / (Zeff + Z1)
    # Again correction needed, since this is the transm coeff for the VELOCITY POTENTIAL (see pdf)
    # We correct to obtain the displacement / velocity coefficient:
    Tsf = Tsf * c_long_liquid / c_shear_solid
    # now Tdf is the displacement( or velocity) transmission coeff, like in Auld, for shear wabes
    VTransm_shear = Tsf

    # now computer the power transmissie coeff for shear(using impedances and angles):
    Tspow = np.abs(Tsf) ** 2 * (rho_solid * c_shear_solid) / (rho_liquid * c_long_liquid) * cosbeta_shear / np.cos(alpha_in)
    # Note: Tspow is real below the critical angle.Above the critical angle is cosbeta_shear zuiver pure imaginairy,
    # hence also Tspow is pure imaginairy.In that case there is no power transmission(evanescent  waves), so
    # we can take the real part of Tspow to describe the actual power transmission:
    Transm_shear = np.real(Tspow)

    # phase of reflection coeff:
    Ph_refl = np.angle(VRefl)
    # phase of transmission coeff longit
    Ph_tr_long = np.angle(VTransm_long)
    # phase of transmission coeff shear
    # Ph_tr_shear = 180 * np.angle(VTransm_shear) / np.pi;
    Ph_tr_shear = np.angle(VTransm_shear)

    return Refl, Transm_long, Transm_shear, Ph_refl, Ph_tr_long, Ph_tr_shear


def refraction(v_in,n,c_in,c_out):
    # refraction with incoming direction v_in, normal n, velocities c_in and c_out v_in, n and v_out
    # must be normalised to length poldir gives the VERTICAL polarization direction(only
    # useful for SHEAR waves)
    v_in = v_in / LA.norm(v_in)
    n = n / LA.norm(n)
    if np.dot(v_in, n) < 0:
       n = -n # n must be directed into out
    nbr = c_out / c_in
    a = v_in - (np.dot(v_in, n))*n
    na = LA.norm(a)
    if na == 0: # incoming ray parallel to normal, i.e.perpendicular to surface
       refr = True
       v_out = v_in
    # in this case the polarization direction of the reflected shear wave is undefined,
    # but the amplitude of this reflected wave is zero
       poldir = np.array([0, 0, 0]) # arbitrary value
    else: # incoming ray not perpendicular to surface
         sintheta2 = nbr * na; # Snellius: sintheta1 / sintheta2 = c_in / c_out
    if sintheta2 <= 1:
        #refraction
       costheta2 = np.sqrt(1 - sintheta2 ** 2)
       refr = True
       v_out = costheta2 * n + nbr * a
       poldir = -sintheta2 * n + costheta2 / na * a
    else: # no refraction
        refr = False
        v_out = np.array([0, 0, 0])
        poldir = np.array([0, 0, 0]) # arbitrary value

    return refr, v_out, poldir


def B2MCoef_long(alpha_in,rho_liquid,rho_solid,c_long_liquid,c_long_solid,c_shear_solid):
    # BEREKENING REFLECTIE TRANSMISSIE COEFF VAN SOLID(LONGIT INPUT, LONGIT en SHEAR REFLECTIE) NAAR LIQUID(LONGIT
    # OUTPUT) MET MATRIX VERGEL OPLOSSEN VOLGENS AULD % alpha_in in radialen !!

    # Lamé pars solid:
    mu_solid = rho_solid * c_shear_solid ** 2
    lambda_solid = rho_solid * c_long_solid ** 2 - 2 * mu_solid

    # Lamé pars liquid:
    lambda_liquid = rho_liquid * c_long_liquid ** 2

    # wave numbers: k c = omega = 2pif, dus k = 2pif / c f = 10 ^ 4; % % f = 1.4MHz value of f
    # arbitrary, select f value such that elements in N matrix of approx same size:
    f = 1e-4;
    k_long_liquid = 2 * np.pi * f / c_long_liquid
    k_long_solid = 2 * np.pi * f / c_long_solid
    k_shear_solid = 2 * np.pi * f / c_shear_solid

    # compute beta_long: Snellius: sin(alpha_in) / sin(beta_long) = c_long_solid / c_long_liquid;
    sinbeta_long = np.sin(alpha_in) * c_long_liquid / c_long_solid
    cosbeta_long = cmath.sqrt(1 - sinbeta_long ** 2); # pure imag if alpha_in > alpha_long_crit

    # compute beta_shear: Snellius: sin(alpha_in) / sin(beta_shear) = c_long_solid / c_shear_solid;
    sinbeta_shear = np.sin(alpha_in) * c_shear_solid / c_long_solid
    cosbeta_shear = cmath.sqrt(1 - sinbeta_shear ** 2); #pure imag if alpha_in > alpha_shear_crit


    # refraction leads to longit wave, reflection to longit and shear wave

    # matrix:
    N = np.matrix([[np.cos(alpha_in), sinbeta_shear, cosbeta_long],
                       [(lambda_solid + 2 * mu_solid * (np.cos(alpha_in)) ** 2) * k_long_solid, 2 * mu_solid * k_shear_solid * cosbeta_shear * sinbeta_shear,
                        -lambda_liquid * k_long_liquid],
                       [2 * mu_solid * k_long_solid * np.cos(alpha_in) * np.sin(alpha_in),
                         mu_solid * k_shear_solid * (sinbeta_shear ** 2 - cosbeta_shear ** 2), 0]])

    # right hand side:
    b = np.matrix([[np.cos(alpha_in)],
                  [-k_long_solid * (lambda_solid + 2 * mu_solid * (np.cos(alpha_in)) ** 2)],
                  [2 * mu_solid * k_long_solid * np.cos(alpha_in) * np.sin(alpha_in)]])

    # oplossen amplitude coeff:
    xx = N ** (-1) * b
    # power coeff: coeff zijn altijd real, er is geen critical angle
    Refl_long = np.abs(xx[0][0]) ** 2
    Refl_shear = np.abs(xx[1][0]) ** 2 * c_shear_solid / c_long_solid * cosbeta_shear / np.cos(alpha_in)
    Transm_long = np.abs(xx[2][0]) ** 2 * c_long_liquid * rho_liquid / (c_long_solid * rho_solid) * cosbeta_long / np.cos(alpha_in)

    return Refl_long, Refl_shear, Transm_long


def reflection2(v_in, n, c1, c2):
    # reflection of wave with incoming direction v_in, normal n, speed c1 outgoing
    # wave with speed c2, pol_direction gives vertical polarization direction in case reflected wave is shear
    # v_in, n and v_out must be normalised to length
    v_in = v_in / LA.norm(v_in)
    n = n / LA.norm(n)
    if np.dot(v_in, n)>0:
        n = -n # n must be directed towards incoming and reflected wave
    nbr = c2 / c1
    a = v_in - (np.dot(v_in, n))*n
    sinalpha_in = LA.norm(a)
    sinalpha_out = nbr * sinalpha_in
    if sinalpha_in == 0:
        # incoming ray parallel to normal, i
        reflect_possible = True
        v_out = v_in
        # in this case the polarization direction of the reflected(shear) wave is undefined, but
        # the amplitude of this reflected wave is zero
        pol_direction = np.array([0, 0, 0])
    else: # incoming ray not perpendicular to surface
        if sinalpha_out > 1:
            reflect_possible = False
            v_out = np.array([0, 0, 0])
            pol_direction = np.array([0, 0, 0])
        else:
            reflect_possible = True;
            cosalpha_out = cmath.sqrt(1 - sinalpha_out ** 2)
            v_out = cosalpha_out * n + nbr * a
            pol_direction = -sinalpha_out * n + cosalpha_out / sinalpha_in * a;

    return reflect_possible,v_out, pol_direction


def B2MCoef_shear(alpha_in,rho_liquid,rho_solid,c_long_liquid,c_long_solid,c_shear_solid):
    # BEREKENINGREFLECTIE TRANSMISSIE COEFF VAN SOLID(SHEAR INPUT VERTICALPOLARIZATION, LONGIT
    # en SHEAR REFLECTIE) NAAR LIQUID(LONGIT OUTPUT) MET MATRIX VERGEL OPLOSSEN VOLGENS AULD
    # LEVERT AMPLITUDE REFL / TRANSMISSIE COEFF(PLOT 1) EN POWER REFL EM TRANSMISSIE COEFF(
    # versie volgens Auld Voorbeeld: solid naar liquid / Marrow of Aluminium naar Water
    # alpha_in in radialen !!

    # Lamé pars solid:
    mu_solid = rho_solid * (c_shear_solid ** 2)
    lambda_solid = rho_solid * (c_long_solid ** 2) - 2 * mu_solid

    # Lamé pars liquid:
    lambda_liquid = rho_liquid * (c_long_liquid ** 2)

    # wave numbers: k c = omega = 2 pi f, dus k = 2 pi f / c
    # f = 10 ^ 4; % % f = 1.4 MHz
    # value of f arbitrary, select f value such that elements in N matrix of approx same size:
    f = 1e-4;
    k_long_liquid = 2 * np.pi * f / c_long_liquid
    k_long_solid = 2 * np.pi * f / c_long_solid
    k_shear_solid = 2 * np.pi * f / c_shear_solid

    # computecritical incoming angles:  alpha_long1_crit = asin(c_shear_solid / c_long_solid);
    # kritieke hoek voor longit reflectie naar solid terug
    # alpha_long1_crit_degr = alpha_long1_crit * 180 / pi;

    # alpha_long2_crit = asin(c_shear_solid / c_long_liquid);
    # c_shear_solid / c_long_liquid > 1
    # dus beta_l(long transmissie) < alpha_in(shear incoming) dus geen kritieke hoek voor longit
    # transmissie naar liquid
    # alpha_long2_crit_degr = alpha_long_crit * 180 / pi;

    # compute alpha_long(reflectie longit): Snellius: sin(alpha_in) / sin(alpha_long) = c_shear_solid / c_long_solid;
    sinalpha_long = np.sin(alpha_in) * c_long_solid / c_shear_solid
    cosalpha_long = cmath.sqrt(1 - sinalpha_long ** 2) # pure imag if alpha_in > alpha_long1_crit

    # compute beta_long(refractie longit): Snellius: sin(alpha_in) / sin(beta_long) = c_shear_solid / c_long_liquid
    sinbeta_long = np.sin(alpha_in) * c_long_liquid / c_shear_solid
    cosbeta_long = np.sqrt(1 - sinbeta_long ** 2); # always real

    # refraction leads to longit wave, reflection to longit and shear(vertically polarized) wave
    # matrix:
    N = np.matrix([[np.sin(alpha_in), cosalpha_long, cosbeta_long],
                   [2 * mu_solid * k_shear_solid * np.sin(alpha_in) * np.cos(alpha_in),
                    2 * mu_solid * k_long_solid * cosalpha_long ** 2 + lambda_solid * k_long_solid,
                    -lambda_liquid * k_long_liquid],
                   [k_shear_solid * ((np.sin(alpha_in)) ** 2 - (np.cos(alpha_in)) ** 2),
                    2 * k_long_solid * sinalpha_long * cosalpha_long, 0]])

    # right hand side:
    b = np.matrix([[-np.sin(alpha_in)],
                   [2 * mu_solid * k_shear_solid * np.sin(alpha_in) * np.cos(alpha_in)],
                   [-k_shear_solid * ((np.sin(alpha_in)) ** 2 - (np.cos(alpha_in)) ** 2)]])

    # oplossen amplitude coeff:
    xx = N ** (-1) * b
    # power coeff:
    Refl_long = np.real(np.abs(xx[0][0]) ** 2)
    Refl_shear = np.real(np.abs(xx[1][0]) ** 2 * c_long_solid  / c_shear_solid * cosalpha_long / np.cos(alpha_in))
    Transm_long = np.real(np.abs(xx[2][0]) ** 2 * c_long_liquid * rho_liquid / (
            c_shear_solid * rho_solid) * cosbeta_long / np.cos(alpha_in))

    return Refl_long,Refl_shear,Transm_long


def reflection3(v_in,n):
    # reflection of shear wave with incoming direction v_in, polarization direction pol, normal n,
    # computes direction and pol directions of the incoming and reflected outgoing shear wave
    # v_in, n and v_out must be normalised to length
    v_in = v_in / LA.norm(v_in)
    n = n / LA.norm(n)
    if np.dot(v_in, n) > 0:
        n = -n # n must be directed towards incoming and reflected wave
    a = v_in - (np.dot(v_in, n))*n # hor component of n
    na = LA.norm(a)
    if na == 0: # incoming ray parallel to normal
        v_out = v_in
        # in this case the polarization direction of the reflected shear wave is undefined,
        # but the amplitude of this reflected wave is zero
        vert_pol_dir = np.array([0, 0, 0]) # arbitrary value
        hor_pol_dir = np.array([0, 0, 0]) # arbitrary value
    else: # incomingraynot perpendicular to surface
        sinalpha_in = na
        cosalpha_in = cmath.sqrt(1 - sinalpha_in ** 2)
        v_out = -v_in + 2 * a
        vpol_dir_in = sinalpha_in * n + cosalpha_in / na * a
        vpol_dir_refl = -sinalpha_in * n + cosalpha_in / na * a
        hor_pol_dir = np.cross(a, n) / sinalpha_in
    return v_out, vpol_dir_in,vpol_dir_refl, hor_pol_dir

