# -*- coding: utf-8 -*-
from Vertice import Vertice
import itertools


class Edge():

    def __init__(self, v1=Vertice(), v2=Vertice()):
        self.vertices = (v1, v2)

    def __eq__(self, other):
        """
        The two edges (u, v) and (v, u) are the same.
        In other words, the pair is not ordered.
        :param other: Edge
        :return: bool
        """
        return set(self.vertices) == set(other.vertices)

    def is_parallel(self, other):
        """
        Edges that have the same end vertices are parallel.
        :param other: Edge
        :return: bool
        """
        return self == other

    def is_loop(self):
        """
        An edge of the form (v, v) is a loop.
        :return: bool
        """
        if len(self.vertices) == 1:
            return True
        return self.vertices[0] == self.vertices[1]

    def is_adjacent(self, other):
        """
        Edges are adjacent if they share a common end vertex.
        :param other: Edge
        :return: bool
        """
        return any([x[0] == x[1] for x in itertools.product(self.vertices, other.vertices)])

    def is_pendant(self):
        """
        An edge that has a pendant vertex as an end vertex is a pendant edge
        :return: bool
        """
        return any([vertice.is_pendant() for vertice in self.vertices])

    def __hash__(self):
        return hash(self.vertices)
