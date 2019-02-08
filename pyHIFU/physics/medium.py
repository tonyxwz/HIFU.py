import numpy as np
from pyHIFU.io.config import readjson
from pyHIFU.geometric.surfaces import Plane
from pyHIFU.geometric.volumes import Ball, Cuboid, Cylinder
from pyHIFU.physics import LIQUID, SOLID, SHEAR, LONGITUDINAL

SHAPE_FUNC_DICT = {'ball': Ball, 'cuboid': Cuboid, 'cylinder': Cylinder}

markoil_properties = {

}

class Material(object):
    """ only physics properties here """
    def __init__(self, material_name=None, state=SOLID, density=None,
            cL=None, cS=None,
            attenuationL=None, attenuationS=None,
            absorption=None, thermal_conductivity=None, heat_capacity=None, **kw):
        self.material_name = material_name
        self.state = state
        self.density = density
        c = [cL]  # velocity
        attenuation = [attenuationL]
        if self.state == SOLID:
            # SOLID == 1
            # self.c[SHEAR], self.c[LONGITUDINAL]
            c.append(cS)
            attenuation.append(attenuationS)
        self.c = np.array(c)
        self.attenuation = np.array(attenuation)
        self.absorption = absorption
        self.thermal_conductivity = thermal_conductivity  #k
        self.heat_capacity = heat_capacity  #cp


    @property
    def Z(self):  # impedence z = c * rho
        d = self.c * self.density
        return d

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

        C = omega**2*rho/(alpha**2+ k**2)
        D = np.sqrt(2)*C/(speed * np.sqrt(rho))
        p_1 = D**2-C
        p_2 = np.sqrt(C**2-p_1**2)/omega
        return p_1, p_2


class Medium(Material):
    def __init__(self, material_name=None, state=None, density=None,
            cL=None, cS=None,
            attenuationL=None, attenuationS=None,
            absorption=None, thermal_conductivity=None, heat_capacity=None,
            med_name=None, med_idx=None,
            geoinfo=None, is_init_med=False, **kw):

        super().__init__(material_name=material_name, state=state, density=density,
                         cL=cL, cS=cS,
                         attenuationL=attenuationL, attenuationS=attenuationS,
                         absorption=absorption, thermal_conductivity=thermal_conductivity, heat_capacity=heat_capacity)
        self.idx = med_idx
        self.name = med_name
        if is_init_med:
            """ give the boundary """
            pass
        else:
            shape_type = geoinfo['shape_type']
            self.shape = SHAPE_FUNC_DICT[shape_type](**geoinfo['parameters'])


    @staticmethod
    def create_markoil(markoil_parameter):
        return Medium(**markoil_parameter)

    @staticmethod
    def create_muscle(muscle_parameter):
        return Medium(**muscle_parameter)

    @staticmethod
    def create_bone(bone_parameter):
        return Medium(**bone_parameter)

    @staticmethod
    def create_liquid(liquid_parameter):
        return Medium(**liquid_parameter)
    @staticmethod
    def create_solid(solid_parameter):
        return Medium(**solid_parameter)


class InitMedium(Material):
    """ Initial medium e.g. markoil and lossless
    """
    def __init__(self, material_name='', state=LIQUID, density=0,
            cL=0, cS=0,
            attenuationL=0, attenuationS=0,
            absorption=0, thermal_conductivity=0, heat_capacity=0,
            med_name=None, med_idx=0,
            boundary=None, **kw):

        super().__init__(material_name=material_name, state=state, density=density,
                         cL=cL, cS=cS,
                         attenuationL=attenuationL, attenuationS=attenuationS,
                         absorption=absorption, thermal_conductivity=thermal_conductivity, heat_capacity=heat_capacity)
        self.boundary = Plane(**boundary)
        self.is_initial = True
        self.shape = [self.boundary]  # in order to be compatible with Medium
        self.id = -1

    @staticmethod
    def from_config(config):
        pass


class MediaComplex(list):
    """
    a list of all the media that is present in the HIFU system
    use adjacency matrix to find neighbours for acoustics related calculation
    """

    # def __init__(self, media_list, adj_mtx):
    def __init__(self, config_json=None, config_file_path=None):
        super().__init__()
        if config_json is None:
            config_json = readjson(json_path=config_file_path)
        # self.__construct_init_medium(config_json["init_medium"])
        self.__construct_other_media(config_json["medium_list"])
        # self.adj_mtx = np.eye(len(self))

    @staticmethod
    def from_config(config):
        pass

    def __construct_init_medium(self, init_medium_config):
        material_name = "markoil"
        state = LIQUID
        density = 1070
        cL = 1430
        attenuationL = 1.04
        absorption=1.04
        thermal_conductivity=0.5
        heat_capacity=4200

        # config should contain the  information of an infinite plane
        init_med = InitMedium(material_name=material_name, state=state,
                              density=density, cL=cL, attenuationL=attenuationL, absorption=absorption,
                              thermal_conductivity=thermal_conductivity, heat_capacity=heat_capacity,
                              **init_medium_config)
        self.append(init_med)

    def __construct_other_media(self, med_list_json):
        for item in med_list_json:
            phyinfo = item["physics"]
            geoinfo = item["geometry"]
            n = len(self)
            self.append(Medium(**phyinfo, geoinfo=geoinfo, med_idx=n))

        # adj_mtx[i,j] = [1,2] : 1,2 are the indices of the faces of self[i]
        # adjacent to self[j]
        self.adj_mtx = list()

        # now you have the list ready, time to build adjacency matrix
        # traverse all media (double compare loop)
        for i,item1 in enumerate(self):
            l = list()
            for j,item2 in enumerate(self):
                if i == j:
                    l.append([])
                else:
                    l.append(item1.shape.adj_at(item2.shape))
            self.adj_mtx.append(l)

    def find_next(self, med_idx, side_idx):
        """
        find the index of next medium, given the index of current medium
        and outgoing side `side`
        """
        # TODO: if there's no adjacency, it is adjacent to the air.
        r = []
        for i,v in enumerate(self.adj_mtx[med_idx]):
            if side_idx in v:
                r.append(i)
        return r

    def find_next_init(self):
        """ not necessary any more """
        for i,v in enumerate(self):
            for j, f in enumerate(v):
               if self[0].boundary == f:
                   return i,j


if __name__ == "__main__":
    mc = MediaComplex(config_file_path='data/example2.json')
    print(mc.adj_mtx)