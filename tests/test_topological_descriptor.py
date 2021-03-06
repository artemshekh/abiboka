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

    def test_superpendentic_index(self):
        assert_almost_equals(superpendentic_index(self.dimethylpentane), 5.831, places=3)

    def test_petitjean_shape_index(self):
        assert_almost_equals(petitjean_shape_index(self.ethynylbenzoicacid), 0.5, places=3)

    def test_eccentricity(self):
        assert_equals(eccentricity(self.ethynylbenzoicacid), 53)

    def test_average_eccentricity(self):
        assert_almost_equals(average_eccentricity(self.dimethylpentane), 3.429, places=3)

    def test_eccentric(self):
        assert_almost_equals(eccentric(self.ethynylbenzoicacid), 0.777, places=3)

    def test_average_graph_distance_degree(self):
        assert_almost_equals(average_graph_distance_degree(self.ethynylbenzoicacid), 28.364, places=3)

    def test_mean_distance_degree_deviation(self):
        assert_almost_equals(mean_square_distance_index(self.ethynylbenzoicacid), 2.236, places=3)

    def test_unipolarity(self):
        assert_equals(unipolarity(self.ethynylbenzoicacid), 21)

    def test_centralization(self):
        assert_equals(centralization(self.ethynylbenzoicacid), 81)

    def test_variance(self):
        assert_equals(variance(self.dimethylpentane), 6)

    def test_radial_centric_information_index(self):
        assert_almost_equals(radial_centric_information_index(self.ethynylbenzoicacid), 1.868, places=3)

    def test_schultz_topological_index(self):
        assert_almost_equals(schultz_topological_index(self.dimethylpentane), 176, places=3)

    def test_schultz_topological_index_by_valence_degree(self):
        assert_almost_equals(schultz_topological_index_by_valence_degree(self.ethanol), 56, places=3)

    def test_gutman_topological_index(self):
        assert_equals(gutman_topological_index(self.ethanol), 6)

    def test_gutman_topological_index_by_valence_degree(self):
        assert_equals(gutman_topological_index_by_valence_degree(self.ethanol), 22)

    def test_xu_index(self):
        assert_almost_equals(xu_index(self.dimethylpentane), 6.786, places=3)

    def test_mti_index(self):
        assert_equals(mti_index(self.dimethylpentane), 150)

    def test_sh_index_1(self):
        assert_almost_equals(sh_index_1(self.ethanol), -0.470, places=3)

    def test_sh_index_2(self):
        assert_almost_equals(sh_index_2(self.ethanol), 1.897, places=3)

    def test_sh_index_3(self):
        assert_almost_equals(sh_index_3(self.ethanol), -1.685, places=3)

    def test_sh_index_4(self):
        assert_almost_equals(sh_index_4(self.ethanol), 0.106, places=3)

    def test_sh_index_5(self):
        assert_almost_equals(sh_index_5(self.ethanol), -0.927, places=3)

    def test_sh_index_6(self):
        assert_almost_equals(sh_index_6(self.ethanol), 3.951, places=3)

    def test_sh_index_7(self):
        assert_almost_equals(sh_index_7(self.ethanol), 3.775, places=3)

    def test_sh_index_8(self):
        assert_almost_equals(sh_index_8(self.ethanol), 3.638, places=3)

    def test_sh_index_9(self):
        assert_almost_equals(sh_index_9(self.ethanol), 4.718, places=3)

    def test_sh_index(self):
        assert_almost_equals(sh_index(self.ethanol), 1.335, places=3)

    def test_all_path_wiener_index(self):
        assert_equals(all_path_wiener_index(self.ethanol), 4)

    def test_kier_alpha_modified_shape_index_1(self):
        assert_almost_equals(kier_alpha_modified_shape_index_1(self.ethanol), 13.278, places=3)

    def test_kier_alpha_modified_shape_index_2(self):
        assert_almost_equals(kier_alpha_modified_shape_index_2(self.ethanol), 1.868, places=3)

    def test_kier_alpha_modified_shape_index_3(self):
        assert_almost_equals(kier_alpha_modified_shape_index_3(self.dimethylpentane), 6.0, places=3)

    def test_kier_molecular_flexibility_index(self):
        assert_almost_equals(kier_molecular_flexibility_index(self.ethanol), 8.270, places=3)