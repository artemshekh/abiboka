# -*- coding: utf-8 -*-
from calc.matrixes.matrix import Matrix
from collections import defaultdict, Counter
from descriptors.topological import close_shell


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



def path_vector(molecule):
    """
    calculate atomic path number for n
    :param molecule:
    :return:
    """
    molecule = molecule.hydrogen_suppressed
    def dfs(atom):
        used_atom.add(atom)
        bonds = atom.bonds
        for bond in bonds:
            if bond not in bonds_stack:
                for _ in bond:
                    if _ is not atom and _ is not atom_[0]:
                        bonds_stack.append(bond)
                        dct[atom_[0]][len(bonds_stack)] += 1
                        dfs(_)
        if bonds_stack:
            bonds_stack.pop()
    dct = defaultdict(lambda : Counter())
    for atom in molecule.atoms:
        atom_ = [atom]
        used_atom = set()
        bonds_stack = []
        dfs(atom)
    return dct

def path_sequence_matrix(molecule, l=None):
    molecule = molecule.hydrogen_suppressed
    dct = path_vector(molecule)
    if not l:
        l = max([max(v.keys()) for k, v in dct.items()])
    m = Matrix([[0 for y in range(l)] for x in range(len(molecule.atoms))])
    for index, atom in enumerate(molecule.atoms):
        for k, v in dct[atom].iteritems():
            if k <= l:
                m.matrix[index][k-1] = v
    return m

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


close_shell = [1, 3, 11, 19, 37, 55, 87]

def valence_electrones(atom):
    valence_e = 4
    for period,n in enumerate(close_shell):
        if atom.Z < n:
            core_e, valence_e = close_shell[period-1]-1, atom.Z - close_shell[period-1] + 1
            break
    if valence_e:
       return valence_e
    else:
        return atom.Z - 87


def valence_degree(atom):
    # all valence electrons of the ith atom
    vd = valence_electrones(atom)
    for _ in atom.connected_with():
        if atom.Z == 1:
            vd -= 1
    return vd


def valence_degree_(atom):
    # the electronic identity of the atom in terms of both valence electron and core electron counts
    vd = valence_electrones(atom)
    for _ in atom.connected_with():
        if atom.Z == 1:
            vd -= 1
    vd = float(vd)/(atom.Z - valence_e - 1)
    return vd