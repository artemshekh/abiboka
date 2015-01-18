# -*- coding: utf-8 -*-
"""
Autocorrelation index
Broto-Moreau autocorrelation index

"""
import math


from calc.matrixes.matrix import Matrix
from descriptors.matrixes import distance_matrix
from utils.periodic_table import periodic_table


def broto_moreau_index(molecule, k, property):
    descriptor = 0
    molecule = molecule.hydrogen_suppressed
    matrix = molecule.distance_matrix
    for i, row in enumerate(matrix):
        for j, value in enumerate(row):
            if value == k:
                descriptor += periodic_table[molecule.atoms[i].Z][property] *\
                              periodic_table[molecule.atoms[j].Z][property]
    if descriptor:
        return math.log(descriptor/2.0)
    else:
        return 0