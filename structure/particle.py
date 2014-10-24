from utils.constants import ELECTRON_MASS, PROTON_MASS, NEUTRON_MASS, ELECTRON_CHARGE

class Particle(object):
    def __init__(self):
        """
        TODO: common properties of particles
        :return:
        """
        raise NotImplementedError()

    def __eq__(self, other):
        return isinstance(self, other.__class__)


class Proton(Particle):
    def __init__(self, name='Proton', weight=PROTON_MASS, charge=+ELECTRON_CHARGE, spin=0.5):
        self.name = name
        self.weight = weight
        self.charge = charge
        self.spin = spin


class Neutron(Particle):
    def __init__(self, name='Neutron', weight=NEUTRON_MASS, charge=0, spin=0.5):
        self.name = name
        self.weight = weight
        self.charge = charge
        self.spin = spin


class Electron(Particle):
    def __init__(self, name='Electron', weight=ELECTRON_MASS, charge=ELECTRON_CHARGE, spin=0.5):
        self.name = name
        self.weight = weight
        self.charge = charge
        self.spin = spin