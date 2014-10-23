# -*- coding: utf-8 -*-
from collections import Counter

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