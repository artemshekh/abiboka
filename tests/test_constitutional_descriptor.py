import math
from nose.tools import assert_almost_equals, assert_equals


from structure.Atom import Atom
from structure.Bond import Bond
from structure.Molecule import Molecule

from descriptors.constitutional import *

class Constitutional_descriptor_Test():
    def setUp(self):
        m = Molecule()
        m.atoms += [Atom(z=6) for x in range(2)]
        m.atoms += [Atom(z=1) for x in range(6)]
        m.atoms += [Atom(z=8) for x in range(1)]
        self.molecule1 = m



    def test_molecular_weigth(self):
        assert_almost_equals(molecular_weight(self.molecule1), 46.069)

    def test_average_molecluar_weight(self):
        assert_almost_equals(average_molecular_weight(self.molecule1), 5.119, places=3)

    def test_sum_van_der_waals_volume(self):
        assert_almost_equals(sum_van_der_waals_volume(self.molecule1), 4.340, places=3)

    def test_mean_van_der_waals_volume(self):
        assert_almost_equals(mean_van_der_waals_volume(self.molecule1), 0.482, places=3)

    def test_sum_of_atom_polarizability(self):
        assert_almost_equals(sum_of_atom_polarizability(self.molecule1), 5.004, places=3)

    def test_mean_sum_of_atom_polarizability(self):
        assert_almost_equals(mean_sum_atom_polarizability(self.molecule1), 0.556, places=3)

    def test_sum_of_first_ionization_potentials(self):
        assert_almost_equals(sum_of_first_ionization_potentials(self.molecule1), 10.455, places=3)

    def test_mean_sum_of_first_ionization_potentials(self):
        assert_almost_equals(mean_sum_of_first_ionization_potentials(self.molecule1), 1.162, places=3)























