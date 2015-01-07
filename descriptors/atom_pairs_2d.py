# -*- coding: utf-8 -*-
"""
Atom 2d pairs
Topological distance sum between pairs of atoms
Presence pairs of atoms at topological distance
Frequency of atoms at topological distance
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


def frequency_of_atoms_at_topological_distance(molecule, atom1, atom2, distance):
    z1 = periodic_table_by_symbol[atom1]['z']
    z2 = periodic_table_by_symbol[atom2]['z']
    dist_matrix = distance_matrix(molecule)
    molecule = molecule.hydrogen_suppressed
    list1 = [index for index, atom in enumerate(molecule.atoms) if atom.Z == z1]
    list2 = [index for index, atom in enumerate(molecule.atoms) if atom.Z == z2]
    descriptor = 0
    for i in list1:
        for j in list2:
            if i != j and dist_matrix[i][j] == distance:
                descriptor += dist_matrix[i][j]
    s = 0
    for i, row in enumerate(dist_matrix):
        for j, value in enumerate(row):
            if value == 2:
                s += value
    s = s/2
    if z1 == z2:
        descriptor = descriptor/2
    if s:
        descriptor = float(descriptor)/s
    else:
        descriptor = 0
    if presence_at_topological_distance(molecule, atom1, atom2, distance):
        print descriptor
    return descriptor
