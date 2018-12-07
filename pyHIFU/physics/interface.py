class Interface(object):
    def __init__(self, medium1, medium2, normal_vector, surface):
        self.medium1 = medium1
        self.medium2 = medium2
        self.normal_vector = normal_vector
        pass

    def get_ratio(self, ray):
        if self.normal_vector.dot(ray.direction) > 0:
            return medium1