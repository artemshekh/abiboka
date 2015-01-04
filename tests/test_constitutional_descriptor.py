import math
from nose.tools import assert_almost_equals, assert_equals

from fixtures.molecule import ethanol, ethynylbenzoicacid
from descriptors.constitutional import *

class Constitutional_descriptor_Test():

    def setUp(self):
        self.ethanol = ethanol()
        self.ethynylbenzoicacid = ethynylbenzoicacid()

    def test_molecular_weigth(self):
        assert_almost_equals(molecular_weight(self.ethanol), 46.069)

    def test_average_molecluar_weight(self):
        assert_almost_equals(average_molecular_weight(self.ethanol), 5.119, places=3)

    def test_sum_van_der_waals_volume(self):
        assert_almost_equals(sum_van_der_waals_volume(self.ethanol), 4.340, places=3)

    def test_mean_van_der_waals_volume(self):
        assert_almost_equals(mean_van_der_waals_volume(self.ethanol), 0.482, places=3)

    def test_sum_of_atom_polarizability(self):
        assert_almost_equals(sum_of_atom_polarizability(self.ethanol), 5.004, places=3)

    def test_mean_sum_of_atom_polarizability(self):
        assert_almost_equals(mean_sum_atom_polarizability(self.ethanol), 0.556, places=3)

    def test_sum_of_first_ionization_potentials(self):
        assert_almost_equals(sum_of_first_ionization_potentials(self.ethanol), 10.455, places=3)

    def test_mean_sum_of_first_ionization_potentials(self):
        assert_almost_equals(mean_sum_of_first_ionization_potentials(self.ethanol), 1.162, places=3)

    def test_sum_of_sanderson_electronegativity(self):
        assert_almost_equals(sum_of_sanderson_electronegativity(self.ethanol), 8.978, places=3)

    def test_mean_sum_of_sanderson_electronegativity(self):
        assert_almost_equals(mean_sum_of_sanderson_electronegativity(self.ethanol), 0.998, places=3)

    def test_number_of_atoms(self):
        assert_equals(number_of_atoms(self.ethanol), 9)

    def test_number_of_non_hydrogen_atoms(self):
        assert_equals(number_of_non_hydrogen_atoms(self.ethanol), 3)

    def test_number_of_bonds(self):
        assert_equals(number_of_bonds(self.ethanol), 8)

    def test_number_of_non_hydrogen_bonds(self):
        assert_equals(number_of_non_hydrogen_bonds(self.ethanol), 2)

    def test_number_of_single_bonds(self):
        assert_equals(number_of_single_bond(self.ethynylbenzoicacid), 9)

    def test_number_of_double_bonds(self):
        assert_equals(number_of_double_bonds(self.ethynylbenzoicacid), 1)

    def test_number_of_triple_bonds(self):
        assert_equals(number_of_triple_bonds(self.ethynylbenzoicacid), 1)

    def test_number_of_aromatic_bonds(self):
        assert_equals(number_of_aromatic_bonds(self.ethynylbenzoicacid), 6)

    def test_sum_of_conventional_bond_order(self):
        assert_equals(sum_of_conventional_bond_order(self.ethynylbenzoicacid), 17)

    def test_rotational_bond_count(self):
        assert_equals(rotatable_bond_count(self.ethynylbenzoicacid), 1)

    def test_rotational_bond_fraction(self):
        assert_almost_equals(rotatable_bond_fraction(self.ethynylbenzoicacid), 0.059, places=3)






















