from utils.periodic_table import periodic_table
from utils.constants import PROTON_MASS, NEUTRON_MASS, ELECTRON_MASS
"""
First Atom
"""

class Atom(object):
    def __init__(self, z=None, n=None, aromatic=False):
        """
        Initialize
        :param z: number of protons
        :param n: number of neutrons
        :return: self
        """
        self.Z = z
        self.N = n or z
        self.charge = 0
        self.bonds = []
        self.aromatic = aromatic
        try:
            self.properties = periodic_table[z]
        except KeyError:
            print 'Error there is no Atom with such atomic number'
            raise

    def get_mass(self):
        """
        absolute weight of Atom
        :return:
        """
        return self.Z*PROTON_MASS + self.N*NEUTRON_MASS + (self.Z - self.charge) * ELECTRON_MASS

    def get_relative_mass(self):
        """
        Get relative_atomic_weight
        :return:
        """
        return self.properties['relative_atomic_weight']

    def __str__(self):
        return periodic_table[self.Z]['symbol'] + '-' + str(self.N)

    def __repr__(self):
        return self.__str__()

if __name__ == '__main__':
    pass