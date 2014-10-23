# -*- coding: utf-8 -*-
from Vertice import Vertice

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