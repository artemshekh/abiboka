# -*- coding: utf-8 -*-
"""
Walk and path counts
"""

from collections import defaultdict, Counter

from calc.matrixes.matrix import Matrix
from utils.functional import cached
from descriptors.ring_descriptor import cyclomatic_number


@cached
def path_vector(molecule):
    """
    Dont do it on the molecules with large connectivity
    :param molecule:
    :return:
    """
    molecule = molecule.hydrogen_suppressed
    if cyclomatic_number(molecule) > 5:
        return "Big problem molecule"

    def dfs(atom):
        for index, atom_from_stack in enumerate(atom_stack):
            key = [atom_from_stack, atom]
            key.sort()
            key = tuple(key)
            value = atom_stack[index:] + [atom]
            value.sort()
            value = tuple(value)
            path_counter[key].add(value)

        atom_stack.append(atom)
        next_atoms = atom.connected_with()
        for next_atom in next_atoms:
            if next_atom not in atom_stack:
                dfs(next_atom)
        atom_stack.pop()

    visited_atoms = set()
    atom_stack = []
    path_counter = defaultdict(lambda: set())
    for atom in molecule.atoms:
        if atom not in visited_atoms:
            visited_atoms.add(atom)
            dfs(atom)
    final_dct = defaultdict(lambda: Counter())
    for k, path_set in path_counter.items():
        a1, a2 = k
        for path in path_set:
            final_dct[a1][len(path) - 1] += 1
            final_dct[a2][len(path) - 1] += 1
    return final_dct


@cached
def atomic_path_count_vector(molecule, m):
    pvector = path_vector(molecule)
    return [pvector[atom][m] for atom in molecule.hydrogen_suppressed.atoms]


def path_sequence_matrix(molecule, l=None):
    """
    Rows of that matrix - is a vertex path code
    :param molecule:
    :param l:
    :return:
    """
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


@cached
def atomic_path_count_sum(molecule):
    """
    For atom in hydrogen suppresed molecule
    :param molecule:
    :return:
    """
    return [sum(row) for row in path_sequence_matrix(molecule)]


def molecular_path_count(molecule, order=1):
    """
    molecular path order
    :param molecule:
    :param order:
    :return:
    """
    m = path_sequence_matrix(molecule, order)
    return sum([row[-1] for row in m.matrix])/2.0


@cached
def molecular_path_code(molecule):
    m = path_sequence_matrix(molecule)
    print [sum(row) for row in m.transpose()]


@cached
def total_path_count(molecule):
    m = path_sequence_matrix(molecule)
    return len(molecule.atoms) + sum([sum(row) for row in m.matrix])/2.0