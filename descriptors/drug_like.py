# -*- coding: utf-8 -*-
"""
Drug like indexes
"""

def hydrogen_bond_donors(molecule):
    descriptor = 0
    electronegative_atoms = [7, 8]
    for atom in molecule.atoms:
        if atom.Z in electronegative_atoms:
            if any([atom.Z == 1 for atom in atom.connected_with()]):
                descriptor += 1
    return descriptor

def hydrogen_bond_acceptors(molecule):
    descriptor = 0
    acceptors_atom = [7, 8]
    for atom in molecule.atoms:
        if atom.Z in acceptors_atom:
            descriptor += 1
    return descriptor