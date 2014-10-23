import unittest
from calc.graph.Vertice import Vertice
from collections import Counter


class Test_vertice(unittest.TestCase):

    def setUp(self):
        self.isolated_vertice = Vertice()

    def test_initialize(self):
        self.assertEqual(self.isolated_vertice.edges, Counter())
        self.assertTrue(len(self.isolated_vertice.edges)<1)