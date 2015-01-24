# -*- coding: utf-8 -*-
"""
Walk and path counts
"""

from collections import defaultdict, Counter

from calc.matrixes.matrix import Matrix
from utils.functional import cached


@cached
def path_vector(molecule):
    """
    calculate atomic path number for n
    :param molecule:
    :return: dict {Atom: Counter()}
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
    dct = defaultdict(lambda: Counter())
    for atom in molecule.atoms:
        atom_ = [atom]
        used_atom = set()
        bonds_stack = []
        dfs(atom)
    return dct


@cached
def atomic_path_count_vector(molecule, m):
    pvector = path_vector(molecule)
    return [pvector[atom][m] for atom in molecule.hydrogen_suppressed.atoms]





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

def mpc(molecule, order=1):
    """
    molecular path order
    :param molecule:
    :param order:
    :return:
    """
    m = path_sequence_matrix(molecule, order)
    return sum([row[-1] for row in m.matrix])/2.0

def tpc(molecule):
    m = path_sequence_matrix(molecule)
    return len(molecule.atoms) + sum([sum(row) for row in m.matrix])/2.0