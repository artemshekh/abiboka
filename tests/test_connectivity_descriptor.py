# -*- coding: utf-8 -*-
import math
from nose.tools import assert_almost_equals, assert_equals

from fixtures.molecule import ethanol, ethynylbenzoicacid
from descriptors.connectivity import *


class Connectivity_descriptor_Test():

    def setUp(self):
        self.ethanol = ethanol()
        self.ethynylbenzoicacid = ethynylbenzoicacid()

    def test_edge_degree(self):
        assert_equals(edge_degree(self.ethynylbenzoicacid), 28)

    def test_valence_connectivity_index_1(self):
        assert_almost_equals(valence_connectivity_index_1(self.ethanol), 1.023, places=3)