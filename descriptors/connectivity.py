# -*- coding: utf-8 -*-
"""
Connectivity descriptors

Edge degree

"""
from utils.functional import cached

@cached
def edge_degree(molecule):
    return sum([sum(row) for row in molecule.edge_adjacency_matrix])