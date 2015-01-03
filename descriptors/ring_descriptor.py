# -*- coding: utf-8 -*-
"""
Ring descriptors
"""


def cyclomatic_number(molecule):
    molecule = molecule.hydrogen_suppressed
    connectivity = 0
    used = set()
    def dfs(atom):
        used.add(atom)
        for bond in atom.bonds:
            for atom_ in bond:
                if atom_ not in used:
                    dfs(atom_)
    for atom in molecule.atoms:
        if atom not in used:
            dfs(atom)
            connectivity += 1
    c = len(molecule.bonds) - len(molecule.atoms) + connectivity
    return c


