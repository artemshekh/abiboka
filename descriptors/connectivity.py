# -*- coding: utf-8 -*-
"""
Connectivity descriptors

Edge degree
Edge degree count
connectivity index (0-5)

"""
import math
import operator
from collections import Counter

from descriptors.vertex_degree import valence_degree
from utils.functional import cached
from ring_descriptor import cyclomatic_number
from vertex_degree import vertex_degree

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

@cached
def mean_randic_connectivity_index(molecule):
    return connectivity_index_1(molecule)/len(molecule.bonds)

@cached
def average_connectivity_index_0(molecule):
    return connectivity_index_0(molecule)/molecule.size

@cached
def average_connectivity_index_1(molecule):
    return connectivity_index_1(molecule)/molecule.size

@cached
def average_connectivity_index_2(molecule):
    return connectivity_index_2(molecule)/molecule.size

@cached
def average_connectivity_index_3(molecule):
    return connectivity_index_3(molecule)/molecule.size

@cached
def average_connectivity_index_4(molecule):
    return connectivity_index_4(molecule)/molecule.size

@cached
def average_connectivity_index_5(molecule):
    return connectivity_index_5(molecule)/molecule.size

@cached
def valence_connectivity_index_0(molecule):
    return sum([1/math.sqrt(valence_degree(atom)) for atom in molecule.hydrogen_suppressed.atoms])

@cached
def valence_connectivity_index_1(molecule):
    descriptor = 0
    for bond in molecule.hydrogen_suppressed.bonds:
        atoms = [atom for atom in bond]
        descriptor += (1/math.sqrt(valence_degree(atoms[0])))*(1/math.sqrt(valence_degree(atoms[1])))
    return descriptor

@cached
def valence_connectivity_index_2(molecule):
    descriptor = 0
    for bond1 in molecule.hydrogen_suppressed.bonds:
        for bond2 in molecule.hydrogen_suppressed.bonds:
            intersect = bond2 & bond1
            if len(intersect) == 1:
                union = bond1|bond2
                descriptor += 1/math.sqrt(reduce(operator.mul, [valence_degree(atom) for atom in union]))
    return descriptor/2

@cached
def valence_connectivity_index_3(molecule):
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
    return sum([1/math.sqrt(reduce(operator.mul, [valence_degree(atom) for atom in subgraph])) for subgraph in subgraphs])

@cached
def valence_connectivity_index_4(molecule):
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
    return sum([1/math.sqrt(reduce(operator.mul, [valence_degree(atom) for atom in subgraph])) for subgraph in subgraphs])

@cached
def valence_connectivity_index_5(molecule):
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
    return sum([1/math.sqrt(reduce(operator.mul, [valence_degree(atom) for atom in subgraph])) for subgraph in subgraphs])

@cached
def average_valence_connectivity_index_0(molecule):
    return valence_connectivity_index_0(molecule)/molecule.size

@cached
def average_valence_connectivity_index_1(molecule):
    return valence_connectivity_index_1(molecule)/molecule.size

@cached
def average_valence_connectivity_index_2(molecule):
    return valence_connectivity_index_2(molecule)/molecule.size

@cached
def average_valence_connectivity_index_3(molecule):
    return valence_connectivity_index_3(molecule)/molecule.size

@cached
def average_valence_connectivity_index_4(molecule):
    return valence_connectivity_index_4(molecule)/molecule.size

@cached
def average_valence_connectivity_index_5(molecule):
    return valence_connectivity_index_5(molecule)/molecule.size

@cached
def bond_order_connectivity_index_0(molecule):
    return sum([1/math.sqrt(atom.bond_vertex_degree) for atom in molecule.hydrogen_suppressed.atoms])

@cached
def bond_order_connectivity_index_1(molecule):
    descriptor = 0
    for bond in molecule.hydrogen_suppressed.bonds:
        atoms = [atom for atom in bond]
        descriptor += (1/math.sqrt(atoms[0].bond_vertex_degree))*(1/math.sqrt(atoms[1].bond_vertex_degree))
    return descriptor

@cached
def bond_order_connectivity_index_2(molecule):
    descriptor = 0
    for bond1 in molecule.hydrogen_suppressed.bonds:
        for bond2 in molecule.hydrogen_suppressed.bonds:
            intersect = bond2 & bond1
            if len(intersect) == 1:
                union = bond1|bond2
                descriptor += 1/math.sqrt(reduce(operator.mul, [atom.bond_vertex_degree for atom in union]))
    return descriptor/2

@cached
def bond_order_connectivity_index_3(molecule):
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
    return sum([1/math.sqrt(reduce(operator.mul, [atom.bond_vertex_degree for atom in subgraph])) for subgraph in subgraphs])

@cached
def bond_order_connectivity_index_4(molecule):
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
    return sum([1/math.sqrt(reduce(operator.mul, [atom.bond_vertex_degree for atom in subgraph])) for subgraph in subgraphs])

@cached
def bond_order_connectivity_index_5(molecule):
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
    return sum([1/math.sqrt(reduce(operator.mul, [atom.bond_vertex_degree for atom in subgraph])) for subgraph in subgraphs])

@cached
def average_bond_order_connectivity_index_0(molecule):
    return bond_order_connectivity_index_0(molecule)/molecule.size

@cached
def average_bond_order_connectivity_index_1(molecule):
    return bond_order_connectivity_index_1(molecule)/molecule.size

@cached
def average_bond_order_connectivity_index_2(molecule):
    return bond_order_connectivity_index_2(molecule)/molecule.size

@cached
def average_bond_order_connectivity_index_3(molecule):
    return bond_order_connectivity_index_3(molecule)/molecule.size

@cached
def average_bond_order_connectivity_index_4(molecule):
    return bond_order_connectivity_index_4(molecule)/molecule.size

@cached
def average_bond_order_connectivity_index_5(molecule):
    return bond_order_connectivity_index_5(molecule)/molecule.size

@cached
def balaban_distance_connectivity_index(molecule):
    descriptor = 0
    molecule = molecule.hydrogen_suppressed
    distance_matrix = molecule.distance_matrix
    for i, row in enumerate(molecule.adjacency_matrix):
        for j, value in enumerate(row):
            if value == 1:
                descriptor += 1/math.sqrt(sum(distance_matrix[i])*sum(distance_matrix[j]))
    return len(molecule.bonds) * descriptor/(cyclomatic_number(molecule)+1)

@cached
def balaban_distance_connectivity_index_j(molecule):
    balaban_distance_connectivity_index(molecule) * (cyclomatic_number(molecule) + 1)

@cached
def balaban_distance_connectivity_index_g(molecule):
    molecule = molecule.hydrogen_suppressed
    a = molecule.size **2 * (cyclomatic_number(molecule) + 1)/ (molecule.size + cyclomatic_number(molecule) + 1)
    return a *balaban_distance_connectivity_index(molecule)

@cached
def jt_index(molecule):
    descriptor = 0
    molecule = molecule.hydrogen_suppressed
    distance_matrix = molecule.distance_matrix
    for i, row in enumerate(molecule.adjacency_matrix):
        for j, value in enumerate(row):
            if value == 1:
                ti = sum(distance_matrix[i])/ float(vertex_degree(molecule.atoms[i]))
                tj = sum(distance_matrix[j])/ float(vertex_degree(molecule.atoms[j]))
                descriptor += 1/math.sqrt(ti*tj)
    return len(molecule.bonds) * descriptor/(cyclomatic_number(molecule)+1)
