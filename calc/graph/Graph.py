# -*- coding: utf-8 -*-
from collections import Counter
import itertools
from Vertice import Vertice
from Edge import Edge
from utils.matrix import AdjacencyMatrix

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

    def max_degree(self):
        """
        The minimum degree of the vertices in a graph G is denoted δ(G) (= 0 if there is an isolated
        vertex in G). Similarly, we write ∆(G) as the maximum degree of vertices in G.
        :return: int
        """
        return max([vertice.degree() for vertice in self.vertices])

    def min_degree(self):
        """
        self.max_degree
        :return: bool
        """
        return min([vertice.degree() for vertice in self.vertices])

    def is_trivial(self):
        """
        A graph with only one vertex is trivial.
        :return: bool
        """
        return len(self.vertices) == 1

    def is_complex(self):
        """
        A simple graph that contains every possible edge between all the vertices is called a complete graph.
        :return: bool
        """
        return self.is_simple()\
               and\
               all([vertices[0].adjacent_to(vertices[1]) for vertices in itertools.combinations(self.vertices, 2)])

    def is_subgraph(self, other):
        """
        The graph G1 = (V1, E1) is a subgraph of G2 = (V2, E2) if
        1. V1 ⊆ V2 and
        2. Every edge of G1 is also an edge of G2.
        :param other: Graph
        :return: bool
        """
        return self.vertices.issubset(other.vertices) and set(self.edges).issubset(set(other.edges))

    @classmethod
    def create_from_adjacency_matrix(cls, matrix):
        vertices = []
        edges = Counter()
        if not matrix.matrix:
            return cls(vertices, edges)
        else:
            vertices = [Vertice(Counter()) for x in range(matrix.n)]
            print vertices[0] == vertices[1]
            edges = Counter()
            for rownumber, row in enumerate(matrix.matrix):
                for colnumber, column in enumerate(row):
                    if column != '0' and (int(rownumber) <= int(colnumber)):

                        column = int(column)
                        print rownumber, colnumber, column
                        e = Edge(vertices[rownumber], vertices[colnumber])
                        edges[e] += column

                        vertices[rownumber].edges[e] +=column
                        vertices[colnumber].edges[e] +=column
        return cls(vertices, edges)


if __name__ == '__main__':
    filename = '/home/a.shehovtsov/PycharmProjects/abiboka/fixtures/file_matrix'
    g = Graph.create_from_adjacency_matrix(AdjacencyMatrix.from_file(open(filename, 'r')))
    for vertice in g.vertices:
        print vertice, vertice.is_pendant(), vertice.edges
