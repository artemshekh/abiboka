# -*- coding: utf-8 -*-
"""
Molecule class
Internal representation
"""

import operator
import itertools

from calc.graph import Graph
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

if __name__ == '__main__':
    print Molecule()