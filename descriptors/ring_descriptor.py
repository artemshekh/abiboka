# -*- coding: utf-8 -*-
"""
Ring descriptors
"""

def cyclomatic_number(molecule):
    connectivity = 0
    used = set()
    for atom in molecule.atoms:
        if atom not in used:
            dfs(atom, used)
            connectivity += 1
    c = len(molecule.bonds) - len(molecule.atoms) + connectivity
    return c
