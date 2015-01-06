# -*- coding: utf-8 -*-
"""
Atom 2d pairs
"""
import itertools

from utils.periodic_table import periodic_table_by_symbol
from descriptors.matrixes import distance_matrix


def topological_distance_sum(molecule, atom1, atom2):
    z1 = periodic_table_by_symbol[atom1]['z']
    z2 = periodic_table_by_symbol[atom2]['z']
    dist_matrix = distance_matrix(molecule)
    molecule = molecule.hydrogen_suppressed
    list1 = [index for index, atom in enumerate(molecule.atoms) if atom.Z == z1]
    list2 = [index for index, atom in enumerate(molecule.atoms) if atom.Z == z2]
    descriptor = 0

    for i in list1:
        for j in list2:
            if i!=j:
                descriptor += dist_matrix[i][j]
    if z1 == z2:
        descriptor = descriptor/2
    return descriptor

def presence_at_topological_distance(molecule, atom1, atom2, distance):
    z1 = periodic_table_by_symbol[atom1]['z']
    z2 = periodic_table_by_symbol[atom2]['z']
    dist_matrix = distance_matrix(molecule)
    molecule = molecule.hydrogen_suppressed
    list1 = [index for index, atom in enumerate(molecule.atoms) if atom.Z == z1]
    list2 = [index for index, atom in enumerate(molecule.atoms) if atom.Z == z2]
    for i in list1:
        for j in list2:
            if i!=j:
                 if dist_matrix[i][j] == distance:
                     return True
    return False