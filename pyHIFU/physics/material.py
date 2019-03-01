import numpy as np
from cached_property import cached_property

from pyHIFU.physics import LIQUID, LONGITUDINAL, SHEAR, SOLID

PHYSICS_PROPERTIES = {
    "markoil": {
        'material_name': 'markoil',
        'state': LIQUID,
        'density': 1070,
        'cL': 1430,
        'absorption': 1.04,
        'attenuationL': 1.04,
        'heat_capacity': 4200,
        'thermal_conductivity': 0.5
    },
    "lossless": {
        'material_name': 'lossless',
        'state': LIQUID,
        'cL': 1380,
        'density': 1030,
        'attenuationL': 0,
        'absorption': 0
    },
    "muscle": {
        "material_name": "muscle",
        "state": LIQUID,
        'density': 1010,
        'cL': 1537,
        'absorption': None,
        'attenuationL': 5.76,
        'heat_capacity': 3720,
        'thermal_conductivity': 0.537
    },
    'bone': {
        "material_name": "bone",
        "state": SOLID,
        'density': 2025,
        'cL': 3736,
        'cS': 1995,
        'absorption': None,
        'attenuationL': 1.9,
        'attenuationS': 2.8,
        'heat_capacity': 3720,
        'thermal_conductivity': 0.487
    }
}


class Material(object):
    """ only physics properties here """

    def __init__(self, material_name, **kw):
        if len(kw) == 0:
            self.material_name = material_name
            self.state = PHYSICS_PROPERTIES[self.material_name]['state']
            self.density = PHYSICS_PROPERTIES[self.material_name]['density']
            c = [PHYSICS_PROPERTIES[self.material_name]['cL']]  # velocity
            attenuation = [
                PHYSICS_PROPERTIES[self.material_name]['attenuationL']
            ]
            if self.state == SOLID:
                # SOLID == 1
                # self.c[SHEAR], self.c[LONGITUDINAL]
                c.append(PHYSICS_PROPERTIES[self.material_name]['cS'])
                attenuation.append(
                    PHYSICS_PROPERTIES[self.material_name]['attenuationS'])
            self.c = np.array(c)
            self.attenuation = np.array(attenuation)
            self.absorption = PHYSICS_PROPERTIES[
                self.material_name]['absorption']
            self.thermal_conductivity = PHYSICS_PROPERTIES[self.material_name][
                'thermal_conductivity']  #k
            self.heat_capacity = PHYSICS_PROPERTIES[self.material_name][
                'heat_capacity']  #cp

    @cached_property
    def Z(self):  # impedence z = c * rho
        return self.c * self.density

    def FSolvePars(self, ray):
        """
        < Theoretical stuff for reference II >
        Calculates the wave coefficients
        """
        wave_type = ray.wave_type
        speed = self.c[wave_type]
        omega = ray.angular_frequency
        k = omega / speed
        alpha = self.attenuation[wave_type]
        rho = self.density

        C = omega**2 * rho / (alpha**2 + k**2)
        D = np.sqrt(2) * C / (speed * np.sqrt(rho))
        p_1 = D**2 - C
        p_2 = np.sqrt(C**2 - p_1**2) / omega
        return p_1, p_2
