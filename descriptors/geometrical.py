# -*- coding: utf-8 -*-
"""
Geometrical descriptors
Gravitational indexes (rough)

"""

import math
import itertools

from utils.periodic_table import periodic_table


def interatomic_distance(atom1, atom2):
    xx = atom1.coords[0]*atom2.coords[0]
    yy = atom1.coords[1]*atom2.coords[1]
    zz = atom1.coords[2]*atom2.coords[2]
    return math.sqrt(xx + yy + zz)


def gravitational_index_1(molecule, property='relative_atomic_weight'):
    descriptor = 0
    for atom1, atom2 in itertools.combinations(molecule.atoms):
        mm = periodic_table[atom1.Z][property] * periodic_table[atom1.Z][property]
        descriptor += mm/(interatomic_distance(atom1, atom2)**2)
    return descriptor


def gravitational_index_2(molecule, property='relative_atomic_weight'):
    descriptor = 0
    for bond in molecule.bonds:
        atoms = [atom for atom in bond]
        atom1, atom2 = atoms[0], atoms[1]
        mm = periodic_table[atom1.Z][property] * periodic_table[atom1.Z][property]
        descriptor += mm/(interatomic_distance(atom1, atom2)**2)
    return descriptor