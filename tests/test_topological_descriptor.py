# -*- coding: utf-8 -*-
from nose.tools import assert_almost_equals, assert_equals

from fixtures.molecule import ethanol
from descriptors.topological import *


class Topological_descriptor_Test():

    def setUp(self):
        self.ethanol = ethanol()

    def test_first_zafreb_index(self):
        assert_equals(first_zagreb_index(self.ethanol), 6)

    def test_platt_number(self):
        assert_equals(platt_number(self.ethanol), 2)

    def test_connection_number(self):
        assert_equals(connection_number(self.ethanol), 1)

    def test_first_zagreb_index_by_valence_degree(self):
        assert_equals(first_zagreb_index_by_valence_degree(self.ethanol), 30)

    def test_first_zagreb_index_by_kupchik_degree(self):
        assert_almost_equals(first_zagreb_index_by_kupchik_degree(self.ethanol), 38.150, places=3)

    def test_first_zagreb_index_by_madan_degree(self):
        assert_almost_equals(first_zagreb_index_by_madan_degree(self.ethanol), 7.438, places=3)

    def test_second_zagreb_index(self):
        assert_equals(second_zagreb_index(self.ethanol), 4)
