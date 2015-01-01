# -*- coding: utf-8 -*-
"""
Connectivity descriptors

Edge degree

"""
import math
import operator
from collections import Counter

from utils.functional import cached

@cached
def edge_degree(molecule):
    return sum([sum(row) for row in molecule.edge_adjacency_matrix])

@cached
def edge_degree_count(molecule):
    counter = Counter()
    for row in molecule.edge_adjacency_matrix:
        counter[sum(row)] += 1
    return counter

@cached
def connectivity_index_0(molecule):
    return sum([1/math.sqrt(atom.vertex_degree) for atom in molecule.hydrogen_suppressed.atoms])

@cached
def connectivity_index_1(molecule):
    descriptor = 0
    for bond in molecule.hydrogen_suppressed.bonds:
        atoms = [atom for atom in bond]
        descriptor += (1/math.sqrt(atoms[0].vertex_degree))*(1/math.sqrt(atoms[1].vertex_degree))
    return descriptor

@cached
def connectivity_index_2(molecule):
    descriptor = 0
    for bond1 in molecule.hydrogen_suppressed.bonds:
        for bond2 in molecule.hydrogen_suppressed.bonds:
            intersect = bond2 & bond1
            if len(intersect) == 1:
                union = bond1|bond2
                descriptor += 1/math.sqrt(reduce(operator.mul, [atom.vertex_degree for atom in union]))
    return descriptor/2

@cached
def connectivity_index_3(molecule):
    descriptor = 0
    molecule = molecule.hydrogen_suppressed
    def dfs(atom, step):
        used_set.add(atom)
        step += 1
        for bond in atom.bonds:
            if bond not in bond_list:
                atoms = [atom for atom in bond]
                for atom in atoms:
                    if atom not in used_set:
                        bond_list.append(bond)
                        if step == 3:
                            atom_list = []
                            for bond in bond_list:
                                atom_list += [atom for atom in bond]
                            s = list(set(atom_list))
                            s.sort()
                            subgraphs.add(tuple(s))
                            if bond_list:
                                bond_list.pop()
                        elif step < 3:
                            dfs(atom, step)
        if bond_list:
            bond_list.pop()
        step -= 1

    subgraphs = set()
    for atom in molecule.atoms:
        used_set = set()
        bond_list = []
        dfs(atom, 0)
    return sum([1/math.sqrt(reduce(operator.mul, [atom.vertex_degree for atom in subgraph])) for subgraph in subgraphs])

@cached
def connectivity_index_4(molecule):
    descriptor = 0
    molecule = molecule.hydrogen_suppressed
    def dfs(atom, step):
        used_set.add(atom)
        step += 1
        for bond in atom.bonds:
            if bond not in bond_list:
                atoms = [atom for atom in bond]
                for atom in atoms:
                    if atom not in used_set:
                        bond_list.append(bond)
                        if step == 4:
                            atom_list = []
                            for bond in bond_list:
                                atom_list += [atom for atom in bond]
                            s = list(set(atom_list))
                            s.sort()
                            subgraphs.add(tuple(s))
                            if bond_list:
                                bond_list.pop()
                        elif step < 4:
                            dfs(atom, step)
        if bond_list:
            bond_list.pop()
        step -= 1

    subgraphs = set()
    for atom in molecule.atoms:
        used_set = set()
        bond_list = []
        dfs(atom, 0)
    return sum([1/math.sqrt(reduce(operator.mul, [atom.vertex_degree for atom in subgraph])) for subgraph in subgraphs])

@cached
def connectivity_index_5(molecule):
    descriptor = 0
    molecule = molecule.hydrogen_suppressed
    def dfs(atom, step):
        used_set.add(atom)
        step += 1
        for bond in atom.bonds:
            if bond not in bond_list:
                atoms = [atom for atom in bond]
                for atom in atoms:
                    if atom not in used_set:
                        bond_list.append(bond)
                        if step == 5:
                            atom_list = []
                            for bond in bond_list:
                                atom_list += [atom for atom in bond]
                            s = list(set(atom_list))
                            s.sort()
                            subgraphs.add(tuple(s))
                            if bond_list:
                                bond_list.pop()
                        elif step < 5:
                            dfs(atom, step)
        if bond_list:
            bond_list.pop()
        step -= 1

    subgraphs = set()
    for atom in molecule.atoms:
        used_set = set()
        bond_list = []
        dfs(atom, 0)
    return sum([1/math.sqrt(reduce(operator.mul, [atom.vertex_degree for atom in subgraph])) for subgraph in subgraphs])