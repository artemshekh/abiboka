import math
from nose.tools import assert_almost_equals, assert_equals

from fixtures.molecule import ethanol, ethynylbenzoicacid, borane, ammiac, phosphane
from fixtures.molecule import hydrogen_chloride, hydrogen_sulfide, hydrogen_bromide, hydrogen_iodide
from descriptors.constitutional import *

class Constitutional_descriptor_Test():

    def setUp(self):
        self.ethanol = ethanol()
        self.ethynylbenzoicacid = ethynylbenzoicacid()
        self.borane = borane()
        self.ammiac = ammiac()
        self.phosphane = phosphane()
        self.hydrogen_sulfide = hydrogen_sulfide()
        self.hydrogen_chloride = hydrogen_chloride()
        self.hydrogen_bromide = hydrogen_bromide()
        self.hydrogen_iodide = hydrogen_iodide()

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

    def test_number_of_hydrogen_atoms(self):
        assert_equals(number_of_hydrogen_atoms(self.ethanol), 6)

    def test_number_of_boron_atoms(self):
        assert_equals(number_of_boron_atoms(self.borane), 1)

    def test_number_of_carbon_atoms(self):
        assert_equals(number_of_carbon_atoms(self.ethanol), 2)

    def test_number_of_nytrogen_atoms(self):
        assert_equals(number_of_nytrogen_atoms(self.ammiac), 1)

    def test_number_of_oxygen_atoms(self):
        assert_equals(number_of_oxygen_atoms(self.ethanol), 1)

    def test_number_of_phosphorous_atoms(self):
        assert_equals(number_of_phosphorous_atoms(self.phosphane), 1)

    def test_number_of_sulfur_atoms(self):
        assert_equals(number_of_sulfur_atoms(self.hydrogen_sulfide), 1)

    def test_number_of_chlorine_atoms(self):
        assert_equals(number_of_chlorine_atoms(self.hydrogen_chloride), 1)

    def test_number_of_bromine_atoms(self):
        assert_equals(number_of_bromine_atoms(self.hydrogen_bromide), 1)

    def test_number_of_iodine_atoms(self):
        assert_equals(number_of_iodine_atoms(self.hydrogen_iodide), 1)

    def test_number_of_heavy_atoms(self):
        assert_equals(number_of_heavy_atoms(self.ethanol), 3)

    def test_number_of_heteroatoms(self):
        assert_equals(number_of_heteroatoms(self.ethanol), 1)

    def test_number_of_halogen_atoms(self):
        assert_equals(number_of_halogen_atoms(self.hydrogen_chloride), 1)

    def test_percentage_of_hydrogen_atoms(self):
        assert_almost_equals(percentage_of_hydrogen_atoms(self.ethanol), 66.667, places=3)

    def test_percentage_of_carbon_atoms(self):
        assert_almost_equals(percentage_of_carbon_atoms(self.ethanol), 22.222, places=3)

    def test_percentage_of_nytrogen_atoms(self):
        assert_almost_equals(percentage_of_nytrogen_atoms(self.ammiac), 25.0, places=3)

    def test_percentage_of_oxygen_atoms(self):
        assert_almost_equals(percentage_of_oxygen_atoms(self.ethanol), 11.111, places=3)

    def test_percentage_of_halogen_atoms(self):
        assert_almost_equals(percentage_of_halogen(self.hydrogen_iodide), 50.0, places=2)

    def test_number_of_csp3_atoms(self):
        assert_equals(number_of_csp3_carbon_atoms(self.ethanol), 2)

    def test_number_of_csp2_atoms(self):
        assert_equals(number_of_csp2_carbon_atoms(self.ethynylbenzoicacid), 7)

    def test_number_of_csp_atoms(self):
        assert_equals(number_of_csp_carbon_atoms(self.ethynylbenzoicacid), 2)
        