# -*- coding: utf-8 -*-
from nose.tools import assert_almost_equals, assert_equals

from fixtures.molecule import ethanol, ethynylbenzoicacid, dimethylpentane
from descriptors.topological import *


class Topological_descriptor_Test():

    def setUp(self):
        self.ethanol = ethanol()
        self.ethynylbenzoicacid = ethynylbenzoicacid()
        self.dimethylpentane = dimethylpentane()

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

    def test_second_zagreb_index_by_valence_degree(self):
        assert_equals(second_zagreb_index_by_valence_degree(self.ethanol), 24)

    def test_second_zagreb_index_by_kuphik_degree(self):
        assert_almost_equals(second_zagreb_index_by_kupchik_degree(self.ethanol), 44.792, places=3)

    def test_second_zagreb_index_by_madan_degree(self):
        assert_almost_equals(second_zagreb_index_by_madan_degree(self.ethanol), 13.660, places=3)

    def test_overall_modified_zagreb_index_0(self):
        assert_almost_equals(overall_modified_zagreb_index_0(self.ethanol), 2.5, places=3)

    def test_overall_modified_zagreb_ondex_by_valence_degree_0(self):
        assert_almost_equals(overall_modified_zagreb_ondex_by_valence_degree_0(self.ethanol), 1.7, places=3)

    def test_overall_modified_zagreb_index_1(self):
        assert_almost_equals(overall_modified_zagreb_index_1(self.ethanol), 1.0, places=3)

    def test_overall_modified_zagreb_index_by_valence_degree_1(self):
        assert_almost_equals(overall_modified_zagreb_index_by_valence_degree_1(self.ethanol), 0.6, places=3)

    def test_quadratic_index(self):
        assert_almost_equals(quadratic_index(self.ethanol), 0, places=3)

    def test_bertz_branching_index(self):
        assert_almost_equals(bertz_branching_index(self.ethanol), 1, places=3)

    def test_narumi_simple_index(self):
        assert_almost_equals(narumi_simple_index(self.ethanol), 0.693, places=3)

    def test_arithmetic_topological_index(self):
        assert_almost_equals(arithmetic_topological_index(self.ethanol), 1.333, places=3)

    def test_harmonic_narumi_index(self):
        assert_almost_equals(harmonic_narumi_index(self.ethanol), 1.2, places=3)

    def test_garmonic_narumi_index(self):
        assert_almost_equals(geometric_narumi_index(self.ethanol), 0.885, places=3)

    def test_total_structure_connectivity_index(self):
        assert_almost_equals(total_structure_connectivity_index(self.ethanol), 1.201, places=3)

    def test_pogliani_index(self):
        assert_almost_equals(pogliani_index(self.ethanol), 7.0, places=3)

    def test_ramification_index_1(self):
        assert_equals(ramification_index_1(self.ethynylbenzoicacid), 1)

    def test_ramification_index_2(self):
        assert_equals(ramification_index_2(self.dimethylpentane), 1)

    def test_benzene_like_index(self):
        assert_almost_equals(benzene_like_index(self.ethanol), 0.171, places=3)

    def test_polarity_wiener_index(self):
        assert_equals(polarity_wiener_index(self.ethynylbenzoicacid), 13)

    def test_product_of_row_sums(self):
        assert_equals(product_of_row_sums(self.dimethylpentane), 18)

    def test_log_product_of_row_sums(self):
        assert_almost_equals(log_product_of_row_sums(self.dimethylpentane), 2.890, places=3)

    def test_mean_square_distance_index(self):
        assert_almost_equals(mean_square_distance_index(self.dimethylpentane), 1.786, places=3)
