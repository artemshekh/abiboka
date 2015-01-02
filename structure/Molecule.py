# -*- coding: utf-8 -*-
"""
Molecule class
Internal representation
"""

import operator
import itertools

from calc.graph import Graph, path_find
from calc.matrixes.matrix import Matrix
from utils.periodic_table import periodic_table
# TODO! Think about consistency of atom, molecule, and bond class
class Molecule():
    def __init__(self):
        self.atoms = []
        self.bonds = []
        self.descriptor_cache = {}

    _hydrogen_suppressed = None
    _edge_adjacency_matrix = None
    _vertex_degree_matrix = None
    _vertex_zagreb_matrix = None
    _modified_vertex_zagreb_matrix = None
    _additive_adjacency_matrix = None
    _adjacency_matrix = None
    _distance_matrix = None
    _burden_matrix = None
    _multigraph_distance_matrix = None

    def add_atom(self, atom):
        self.atoms.add(atom)

    def add_atoms(self, atoms):
        for atom in atoms:
            self.add_atom(atom)

    def add_bond(self, bond):
        self.bonds.add(bond)

    def add_bonds(self, bonds):
        for bond in bonds:
            self.add_bond(bond)

    def molecular_mass(self):
        return sum([atom.get_relative_mass() for atom in self.atoms])

    def molecular_graph(self):
        return Graph.Graph.from_molecule(self.atoms, self.bonds)

    @property
    def hydrogen_suppressed(self):
        if self._hydrogen_suppressed:
            return self._hydrogen_suppressed
        else:
            atoms, bonds = [], []
        for atom in self.atoms:
            if atom.Z != 1:
                atoms.append(atom)
                bonds_list = []
                for bond in atom.bonds:
                    h_in_bond_ = False
                    for atom_ in bond:
                        if atom_.Z == 1:
                            h_in_bond_ = True
                    if not h_in_bond_:
                        bonds_list.append(bond)
                atom.bonds = bonds_list
        for bond in self.bonds:
            h_in_bond = False
            for atom in bond:
                if atom.Z == 1:
                    h_in_bond = True
            if not h_in_bond:
                bonds.append(bond)
        molecule = Molecule()
        molecule.atoms = atoms
        molecule.bonds = bonds
        self._hydrogen_suppressed = molecule
        return molecule

    @property
    def size(self):
        return len(self.atoms)

    @property
    def edge_adjacency_matrix(self):
        if self._edge_adjacency_matrix:
            return self._edge_adjacency_matrix
        else:
            bonds = self.hydrogen_suppressed.bonds
            n = len(bonds)
            m = [[0 for x in range(n)] for y in range(n)]
            for i, bond1 in enumerate(bonds):
                for j, bond2 in enumerate(bonds):
                    if bond1 is not bond2:
                        if bond1 & bond2:
                            m[i][j] = 1
            self._edge_adjacency_matrix = m
            return self._edge_adjacency_matrix

    @property
    def vertex_degree_matrix(self):
        if self._vertex_degree_matrix:
            return self._vertex_degree_matrix
        else:
            n = self.hydrogen_suppressed.size
            m = [[0 for x in range(n)] for y in range(n)]
            for index,atom in enumerate(self.hydrogen_suppressed.atoms):
                m[index][index] = atom.vertex_degree
            self._vertex_degree_matrix = m
            return self._vertex_degree_matrix

    @property
    def vertex_zagreb_matrix(self):
        if self._vertex_zagreb_matrix:
            return self._vertex_zagreb_matrix
        else:
            n = self.hydrogen_suppressed.size
            m = [[0 for x in range(n)] for y in range(n)]
            for index,atom in enumerate(self.hydrogen_suppressed.atoms):
                m[index][index] = atom.vertex_degree**2
            self._vertex_zagreb_matrix = m
            return self._vertex_zagreb_matrix

    @property
    def modified_vertex_zagreb_matrix(self):
        if self._modified_vertex_zagreb_matrix:
            return self._modified_vertex_zagreb_matrix
        else:
            n = self.hydrogen_suppressed.size
            m = [[0 for x in range(n)] for y in range(n)]
            for index,atom in enumerate(self.hydrogen_suppressed.atoms):
                m[index][index] = 1.0/(atom.vertex_degree**2)
            self._modified_vertex_zagreb_matrix = m
            return self._modified_vertex_zagreb_matrix

    @property
    def adjacency_matrix(self):
        if self._adjacency_matrix:
            return self._adjacency_matrix
        else:
            n = self.hydrogen_suppressed.size
            m = [[0 for x in range(n)] for y in range(n)]
            d = {}
            for index, atom in enumerate(self.hydrogen_suppressed.atoms):
                d[atom] = index
            for bond in self.hydrogen_suppressed.bonds:
                _ = list(bond)
                i, j = d[_[0]], d[_[1]]
                m[i][j], m[j][i] = 1, 1
            self._adjacency_matrix = m
            return self._adjacency_matrix

    @property
    def additive_adjacency_matrix(self):
        if self._additive_adjacency_matrix:
            return self._additive_adjacency_matrix
        else:
            n = self.hydrogen_suppressed.size
            m = [[0 for x in range(n)] for y in range(n)]
            for i, row in enumerate(self.adjacency_matrix):
                for j, value in enumerate(row):
                    if value:
                        m[i][j] = sum(self.adjacency_matrix[j])
            self._additive_adjacency_matrix = m
            return self._additive_adjacency_matrix

    @property
    def distance_matrix(self):
        if self._distance_matrix:
            return self._distance_matrix
        else:
            adjacency_matrix = self.adjacency_matrix
            m = [row[:] for row in adjacency_matrix]
            for i, row in enumerate(m):
                for j, value in enumerate(row):
                    if i!=j and value == 0:
                        m[i][j] = 1000
            for k in range(len(m)):
                for i in range(len(m)):
                    for j in range(len(m)):
                        m[i][j] = min(m[i][j], m[i][k] + m[k][j])
            self._distance_matrix = m
            return self._distance_matrix

    @property
    def burden_matrix(self):
        if self._burden_matrix:
            return self._burden_matrix
        else:
            adjacency_matrix = self.adjacency_matrix
            m = [row[:] for row in adjacency_matrix]
            for i, row in enumerate(m):
                for j, value in enumerate(row):
                    if i==j:
                        m[i][j] = self.atoms[i].Z
                    elif value == 1:
                        for bond in self.bonds:
                            if self.atoms[i] in bond and self.atoms[j] in bond:
                                m[i][j] = bond.conventional_bond_order * 0.1
                    else:
                        m[i][j] = 0.001
            self._burden_matrix = m
            return self._burden_matrix

    @property
    def multigraph_distance_matrix(self):
        if self._multigraph_distance_matrix:
            return self._multigraph_distance_matrix
        else:
            molecule = self.hydrogen_suppressed
            n = self.hydrogen_suppressed.size
            m = [[0 for x in range(n)] for y in range(n)]
            dct = path_find(molecule)
            for i, row in enumerate(m):
                for j, v in enumerate(row):
                    if i!= j:
                        atom1 = molecule.atoms[i]
                        atom2 = molecule.atoms[j]
                        key = [atom1, atom2]
                        key.sort()
                        key = tuple(key)
                        bonds = dct[key]
                        s = 0
                        for bond in bonds:
                            s += 1.0/bond.conventional_bond_order
                        m[i][j], m[j][i] = s, s

            self._multigraph_distance_matrix = m
            return self._multigraph_distance_matrix






if __name__ == '__main__':
    print Molecule()