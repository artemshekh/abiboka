# -*- coding: utf-8 -*-
from nose.tools import assert_almost_equals, assert_equals

from fixtures.molecule import ethanol, ethynylbenzoicacid, dimethylpentane
from descriptors.vertex_degree import *

class Vertex_descriptor_Test():

    def setUp(self):
        self.ethanol = ethanol()
        self.ethynylbenzoicacid = ethynylbenzoicacid()
        self.dimethylpentane = dimethylpentane()

    def test_vertex_degree(self):
        assert_equals(vertex_degree(self.ethanol.atoms[1]), 2)

    def test_extended_connectivity(self):
        assert_equals(extended_connectivity(self.dimethylpentane.atoms[1]), 4)

    def test_dual_degree(self):
        assert_almost_equals(dual_degree(self.dimethylpentane.atoms[1]), 1.333, places=3)

    def test_number_of_sigma_electrones(self):
        assert_equals(number_of_sigma_electrones(self.dimethylpentane.atoms[1]), 4)

    def test_number_of_bonded_hydrogen(self):
        assert_equals(number_of_bonded_hydrogen(self.dimethylpentane.atoms[1]), 1)

    def test_valence_degree(self):
        assert_equals(valence_degree(self.ethanol.atoms[-1]), 5)

    def test_bond_vertex_degree(self):
        assert_equals(bond_vertex_degree(self.ethanol.atoms[0]), 1)

    def test_atomic_multigraph_factor(self):
        assert_almost_equals(atomic_multigraph_factor(self.ethynylbenzoicacid.atoms[3]), 1, places=1)
