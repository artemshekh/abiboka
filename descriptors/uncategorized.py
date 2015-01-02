# -*- coding: utf-8 -*-
"""
Can't classify category
"""
from utils.functional import cached

@cached
def total_adjacency_index(molecule):
    return 2 * len(molecule.bonds)

@cached
def average_total_adjacency_index(molecule):
    return float(total_adjacency_index(molecule))/molecule.size

@cached
def density_index(molecule):
    return average_total_adjacency_index(molecule)/molecule.size

