# -*- coding: utf-8 -*-
from collections import Counter
from itertools import product


class Graph(object):

    def __init__(self, vertices=set(), edges=Counter()):
        """
        Graph is a pair of sets (V, E), where V is the set of vertices and E is the set of
        edges, formed by pairs of vertices. E is a multiset, in other words, its elements can occur more
        than once so that every element has a multiplicity
        :param vertices:
        :param edges:
        :return:
        """
        self.vertices = vertices
        self.edges = edges

    def validate(self):
        pass
        return True

    def is_simple(self):
        """
        A graph is simple if it has no parallel edges or loops
        :return: bool
        """
        if any([edge.is_loop() for edge in self.edges]):
            return False
        if self.edges.most_common(1) > 1:
            return False
        return True

    def is_empty(self):
        """
        A graph with no edges (i.e. E is empty) is empty.
        :return: bool
        """
        return len(self.edges.elements()) == 0

    def is_null(self):
        """
         A graph with no vertices (i.e. V and E are empty) is a null graph.
        :return: bool
        """
        return self.is_empty() and not self.vertices

    def is_trivial(self):
        """
        A graph with only one vertex is trivial.
        :return: bool
        """
        return len(self.vertices) == 1