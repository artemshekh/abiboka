# -*- coding: utf-8 -*-
"""
All type of matrices which can be helpful for calculation of descriptor
Adjacency matrix
Edge adjacency matrix
Vertex degree matrix
Atom connectivity matrix
Augmented adjacency matrix
Degree adjacency matrix
Distance sum connectivity matrix
Zagreb matrix
Edge zagreb matrix
Modified edge zagreb matrix

"""
import math
import operator

from calc.matrixes.matrix import Matrix
from utils.functional import cached
from utils.periodic_table import periodic_table
from descriptors.vertex_degree import vertex_degree, edge_degree


@cached
def adjacency_matrix(molecule):
    molecule = molecule.hydrogen_suppressed
    size = molecule.size
    matrix = [[0 for x in range(size)] for y in range(size)]
    d = {}
    for index, atom in enumerate(molecule.atoms):
        d[atom] = index
    for bond in molecule.bonds:
        atoms = list(bond)
        i, j = d[atoms[0]], d[atoms[1]]
        matrix[i][j], matrix[j][i] = 1, 1
    return Matrix(matrix)


@cached
def edge_adjacency_matrix(molecule):
    molecule = molecule.hydrogen_suppressed
    bonds = molecule.bonds
    n = len(bonds)
    matrix = [[0 for x in range(n)] for y in range(n)]
    for i, bond1 in enumerate(bonds):
        for j, bond2 in enumerate(bonds):
            if bond1 is not bond2:
                if bond1 & bond2:
                    matrix[i][j] = 1
    return Matrix(matrix)


def generalized_vertex_degree_matrix(molecule, func):
    size = molecule.hydrogen_suppressed.size
    matrix = [[0 for x in range(size)] for y in range(size)]
    for i, atom in enumerate(molecule.hydrogen_suppressed.atoms):
        matrix[i][i] = func(vertex_degree(atom))
    return Matrix(matrix)


@cached
def vertex_degree_matrix(molecule):
    return generalized_vertex_degree_matrix(molecule, lambda x: x)


@cached
def vertex_zagreb_matrix(molecule):
    return generalized_vertex_degree_matrix(molecule, lambda x: x**2)


@cached
def modified_vertex_zagreb_matrix(molecule):
    return generalized_vertex_degree_matrix(molecule, lambda x: 1.0/(x**2))

@cached
def atom_connectivity_matrix(molecule):
    """
    Atom connectivity matrix C
    :param molecule:
    :return:
    """
    molecule = molecule.hydrogen_suppressed
    size = molecule.size
    matrix = [[0 for x in range(size)] for y in range(size)]
    d = {}
    for index, atom in enumerate(molecule.atoms):
        d[atom] = index
    for bond in molecule.bonds:
        atoms = list(bond)
        i, j = d[atoms[0]], d[atoms[1]]
        matrix[i][j], matrix[j][i] = bond.conventional_bond_order, bond.conventional_bond_order
    return Matrix(matrix)

@cached
def augmented_adjacency_matrix(molecule, property):
    molecule = molecule.hydrogen_suppressed
    size = molecule.size
    matrix = [[0 for x in range(size)] for y in range(size)]
    d = {}
    for index, atom in enumerate(molecule.atoms):
        d[atom] = index
    for bond in molecule.bonds:
        atoms = list(bond)
        i, j = d[atoms[0]], d[atoms[1]]
        matrix[i][j], matrix[j][i] = 1, 1
    for i, atom in enumerate(molecule.atoms):
        matrix[i][i] = periodic_table[atom.Z][property]
    return Matrix(matrix)

@cached
def degree_adjacency_matrix(molecule):
    molecule = molecule.hydrogen_suppressed
    size = molecule.size
    matrix = [[0 for x in range(size)] for y in range(size)]
    d = {}
    for index, atom in enumerate(molecule.atoms):
        d[atom] = index
    for bond in molecule.bonds:
        atoms = list(bond)
        i, j = d[atoms[0]], d[atoms[1]]
        value = 1.0/math.sqrt(vertex_degree(atoms[0])*vertex_degree(atoms[1]))
        matrix[i][j], matrix[j][i] = value, value
    return Matrix(matrix)

@cached
def distance_sum_connectivity_matrix(molecule):
    adj_matrix = adjacency_matrix(molecule)
    dist_matrix = distance_matrix(molecule)
    matrix = [row[:] for row in adj_matrix.matrix]
    for i, row in enumerate(adj_matrix.matrix):
        for j, value in enumerate(row):
            if value:
                v = 1.0/sum(dist_matrix.matrix[i])*sum(dist_matrix.matrix[j])
                matrix[i][j] = v
    return Matrix(matrix)

def zagreb_matrix(molecule, exp):
    molecule = molecule.hydrogen_suppressed
    size = molecule.size
    matrix = [[0 for x in range(size)] for y in range(size)]
    d = {}
    for index, atom in enumerate(molecule.atoms):
        d[atom] = index
    for bond in molecule.bonds:
        atoms = list(bond)
        i, j = d[atoms[0]], d[atoms[1]]
        value = math.pow(vertex_degree(atoms[0])*vertex_degree(atoms[1]), exp)
        matrix[i][j], matrix[j][i] = value, value
    return Matrix(matrix)

@cached
def edge_zagreb_matrix(molecule):
    return zagreb_matrix(molecule, 1)

@cached
def modified_edge_zagreb_matrix(molecule):
    return zagreb_matrix(molecule, -1)

def extended_vertex_adjacency_matrix(molecule):
    molecule = molecule.hydrogen_suppressed
    size = molecule.size
    matrix = [[0 for x in range(size)] for y in range(size)]
    d = {}
    for index, atom in enumerate(molecule.atoms):
        d[atom] = index
    for bond in molecule.bonds:
        atoms = list(bond)
        i, j = d[atoms[0]], d[atoms[1]]
        v1, v2 = float(vertex_degree(atoms[0])), float(vertex_degree(atoms[1]))
        value = (v1/v2 + v2/v1)/2
        matrix[i][j], matrix[j][i] = value, value
    return Matrix(matrix)

def heteroatom_corrected_extended_adjacency_matrix(molecule):
    molecule = molecule.hydrogen_suppressed
    size = molecule.size
    matrix = [[0 for x in range(size)] for y in range(size)]
    d = {}
    for index, atom in enumerate(molecule.atoms):
        d[atom] = index
    for bond in molecule.bonds:
        atoms = list(bond)
        i, j = d[atoms[0]], d[atoms[1]]
        e1 = periodic_table[atoms[0].Z]['sanderson_electronegativity']
        e2 = periodic_table[atoms[1].Z]['sanderson_electronegativity']
        v1, v2 = float(vertex_degree(atoms[0])), float(vertex_degree(atoms[1]))
        h1, h2 = v1*e1, v2*e2
        value = (h1/h2 + h2/h1)/2
        matrix[i][j], matrix[j][i] = value, value
    for index, atom in enumerate(molecule.atoms):
        matrix[index][index] = periodic_table[atom.Z]['sanderson_electronegativity']
    return Matrix(matrix)

def heteroatom_multiplicity_corrected_extended_adjacency_matrix(molecule):
    molecule = molecule.hydrogen_suppressed
    size = molecule.size
    matrix = [[0 for x in range(size)] for y in range(size)]
    d = {}
    for index, atom in enumerate(molecule.atoms):
        d[atom] = index
    for bond in molecule.bonds:
        atoms = list(bond)
        i, j = d[atoms[0]], d[atoms[1]]
        e1 = periodic_table[atoms[0].Z]['sanderson_electronegativity']
        e2 = periodic_table[atoms[1].Z]['sanderson_electronegativity']
        v1, v2 = float(vertex_degree(atoms[0])), float(vertex_degree(atoms[1]))
        h1, h2 = v1*e1, v2*e2
        rcp = 1.0/bond.conventional_bond_order # reciprocal conventional bond order
        value = ((h1 +1 - rcp)/(h2 + 1 - rcp) + (h2 + 1 - rcp)/(h1 + 1 - rcp))/2
        matrix[i][j], matrix[j][i] = value, value
    for index, atom in enumerate(molecule.atoms):
        matrix[index][index] = periodic_table[atom.Z]['sanderson_electronegativity']
    return Matrix(matrix)

def extended_edge_adjacency_matrix(molecule):
    molecule = molecule.hydrogen_suppressed
    bonds = molecule.hydrogen_suppressed.bonds
    n = len(bonds)
    matrix = [[0 for x in range(n)] for y in range(n)]
    for i, bond1 in enumerate(bonds):
        for j, bond2 in enumerate(bonds):
            if bond1 is not bond2:
                if bond1 & bond2:
                    e1 = edge_degree(bond1)
                    e2 = edge_degree(bond2)
                    value = (e1/e2 + e2/e1)/2
                    matrix[i][j], matrix[j][i] = value, value
    return Matrix(matrix)


@cached
def distance_matrix(molecule):
    adj_matrix = adjacency_matrix(molecule)
    matrix = [row[:] for row in adj_matrix.matrix]
    for i, row in enumerate(matrix):
        for j, value in enumerate(row):
            if i != j and value == 0:
                matrix[i][j] = 1000
    for k in range(len(matrix)):
        for i in range(len(matrix)):
            for j in range(len(matrix)):
                matrix[i][j] = min(matrix[i][j], matrix[i][k] + matrix[k][j])
    return Matrix(matrix)



def molecular_matrix(molecule):
    # molecular geometry
    raise NotImplementedError

def augmented_molecular_matrix(molecule):
    # molecular geometry
    raise NotImplementedError

def z_matrix(molecule):
    # molecular_geometry
    raise NotImlpementedError

def geometry_distance(molecule):
    # molecular_geometry
    raise NotImplementedError

def topographic_matrix(molecule):
    # molecular_geometry
    raise NotImplementedError

def bond_length_weighted_adjacency_matrix(molecule):
    # molecular geometry
    raise NotImplementedError

def neighborhood_geometry_matrix(molecule):
    # molecular geometry
    raise NotImplementedError

def reciprocal_geometry_matrix(molecule):
    # molecular geometry
    raise NotImplementedError

def reciprocal_topographic_matrix(molecule):
    # molecular geometry
    raise NotImplementedError

#geometric distance/topological distance quotient matrix
#topographic distance/topological distance quotient matrix
#topological distance/geometric distance quotient matrix
#topological distance/topographic distance quotient matrix
#distance/distance combined matrices

# WHIM weighted covariance matrices (WWC matrices)
# charge dencity matrix





