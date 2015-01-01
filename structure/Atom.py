from utils.periodic_table import periodic_table
from utils.constants import PROTON_MASS, NEUTRON_MASS, ELECTRON_MASS
"""
First Atom
"""

class Atom(object):
    def __init__(self, z=None, n=None, charge=0, aromatic=False, chiral=False, bracket = False):
        """
        Initialize
        :param z: number of protons
        :param n: number of neutrons
        :return: self
        """
        self.Z = z
        self.N = n or z
        self.charge = charge
        self.chiral = chiral
        self.bonds = []
        self.aromatic = aromatic
        self.bracket = bracket
        try:
            self.properties = periodic_table[z]
        except KeyError:
            print 'Error there is no Atom with such atomic number'
            raise

    @property
    def vertex_degree(self):
        return len(self.bonds)

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

    def connected_with(self):
        _ = []
        for bond in self.bonds:
            for atom in bond:
                if atom is not self:
                    _.append(atom)
        return _

    def __str__(self):
        return periodic_table[self.Z]['symbol'] + '-' + str(self.N)

    def __repr__(self):
        return self.__str__()

if __name__ == '__main__':
    pass