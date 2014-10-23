# -*- coding: utf-8 -*-
from collections import Counter


class Vertice():

    def __init__(self, edges=Counter()):
        self.edges = edges

    def adjacent_to(self, other):
        """
        Two vertices u and v are adjacent if they are connected by an edge, in other words, (u, v)
        is an edge.
        :param other: Vertice
        :return: bool
        """
        return bool(set(self.edges.elements()) & set(other.edges.elements()))

    def degree(self):
        """
        The degree of the vertex v, written as d(v), is the number of edges with v as an end vertex.
        By convention, we count a loop twice and parallel edges contribute separately.
        :return: int
        """
        return len(self.edges)

    def is_pendant(self):
        """
        A pendant vertex is a vertex whose degree is 1
        :return: bool
        """
        print self.degree()
        return self.degree() == 1

    def is_isolated(self):
        """
         An isolated vertex is a vertex whose degree is 0.
        :return: bool
        """
        return self.degree() == 0