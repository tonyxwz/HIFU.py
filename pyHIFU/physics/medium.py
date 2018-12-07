import numpy as np


class Medium(object):
    """ Medium """
    def __init__(self, type=None, name=None,
            speed=None, density=None, attenuation=None, absorption=None,
            frequency=None,
            faces=None, adj_mtx=None, config_file=None):
        self.name = name
        self.type = type
        self.speed = speed
        self.density = density
        self.attenuation = attenuation
        self.absorption = absorption
        self.frequency = frequency
        self.angular_freq = 2 * np.pi * self.frequency

