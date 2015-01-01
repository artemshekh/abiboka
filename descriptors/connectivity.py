# -*- coding: utf-8 -*-
"""
Connectivity descriptors

Edge degree

"""
import math
from collections import Counter

from utils.functional import cached

@cached
def edge_degree(molecule):
    return sum([sum(row) for row in molecule.edge_adjacency_matrix])

@cached
def edge_degree_count(molecule):
    counter = Counter()
    for row in molecule.edge_adjacency_matrix:
        counter[sum(row)] += 1
    return counter

@cached
def connectivity_index_0(molecule):
    return sum([1/math.sqrt(atom.vertex_degree) for atom in molecule.hydrogen_suppressed.atoms])

@cached
def connectivity_index_1(molecule):
    descriptor = 0
    for bond in molecule.hydrogen_suppressed.bonds:
        atoms = [atom for atom in bond]
        descriptor += (1/math.sqrt(atoms[0].vertex_degree))*(1/math.sqrt(atoms[1].vertex_degree))
    return descriptor