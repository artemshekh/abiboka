# -*- coding: utf-8 -*-
from calc.matrixes.matrix import Matrix
from collections import defaultdict, Counter
from descriptors.walk import path_vector


def p3_matrix(molecule):
    matrix = [[0 for y in range(len(molecule.atoms))] for x in range(len(molecule.atoms))]
    d = {}
    for i, _ in enumerate(molecule.atoms):
        d[_] = i
    for i, startatom in enumerate(molecule.atoms):
        walked_atom = set([startatom])
        second_iter_atom = set()
        third_iter_atom = set()
        final_atom = set()
        for bond in startatom.bonds:
            for atom in bond:
                if atom not in walked_atom:
                    second_iter_atom.add(atom)
        for second_atom in second_iter_atom:
            walked_atom.add(second_atom)
            for bond in second_atom.bonds:
                for atom in bond:
                    if atom not in walked_atom:
                        third_iter_atom.add(atom)
        for third_atom in third_iter_atom:
            walked_atom.add(third_atom)
            for bond in third_atom.bonds:
                for atom in bond:
                    if atom not in walked_atom:
                        final_atom.add(atom)
        for atom in final_atom:
            j = d[atom]
            matrix[i][j], matrix[j][i] = 3, 3
    return Matrix(matrix)


def p2_matrix(molecule):
    matrix = [[0 for y in range(len(molecule.atoms))] for x in range(len(molecule.atoms))]
    d = {}
    for i, _ in enumerate(molecule.atoms):
        d[_] = i
    for i, startatom in enumerate(molecule.atoms):
        walked_atom = set([startatom])
        second_iter_atom = set()
        final_atom = set()
        for bond in startatom.bonds:
            for atom in bond:
                if atom not in walked_atom:
                    second_iter_atom.add(atom)
        for second_atom in second_iter_atom:
            walked_atom.add(second_atom)
            for bond in second_atom.bonds:
                for atom in bond:
                    if atom not in walked_atom:
                        final_atom.add(atom)
        for atom in final_atom:
            j = d[atom]
            matrix[i][j], matrix[j][i] = 2, 2
    return Matrix(matrix)







def walk_vector(molecule, order):
    molecule = molecule.hydrogen_suppressed
    dct = defaultdict(lambda : 0)
    def walk(atom, step):
        if step < order:
            step += 1
            for next_atom in atom.connected_with():
                walk(next_atom, step)
        else:
            dct[atom] += 1
    for atom in molecule.atoms:
        walk(atom, 0)
    return dct

