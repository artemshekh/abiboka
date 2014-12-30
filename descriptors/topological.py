# -*- coding: utf-8 -*-
"""
Topological descriptors
"""

def zm1(molecule):
    return sum([len(atom.bonds)**2 for atom in molecule.hydrogen_suppressed().atoms])

def zm1_H(molecule):
    return sum([len(atom.bonds)**2 for atom in moleculeatoms])

