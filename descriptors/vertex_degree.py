# -*- coding: utf-8 -*-
"""
Invariants of atom in molecule

Vertex degree
Extended connectivity
Dual degree
Number of sigma electrones
Number of bonded hydrogens


"""
import math
import operator

from utils.periodic_table import periodic_table


sqrt_2 = math.sqrt(2)
close_shell = [0, 2, 10, 18, 36, 54, 86]


def vertex_degree(atom):
    return sum([atom.Z != 1 for atom in atom.connected_with()])


def extended_connectivity(atom):
    return sum([vertex_degree(atom) for atom in atom.connected_with() if atom.Z != 1])


def dual_degree(atom):
    return float(extended_connectivity(atom))/vertex_degree(atom)


def number_of_sigma_electrones(atom):
    return len([atom for atom in atom.connected_with()])


def number_of_bonded_hydrogen(atom):
    return number_of_sigma_electrones(atom) - vertex_degree(atom)



def valence_electrones(atom):
    return atom.Z - close_shell[periodic_table[atom.Z]['principal_quantum_number']-1]


def valence_degree(atom):
    # all valence electrons of the ith atom
    ve = valence_electrones(atom)
    for atom in atom.connected_with():
        if atom.Z == 1:
            ve -= 1
    return ve


def kier_hall_electronegativity(atom):
    return float(valence_degree(atom) - vertex_degree(atom))/(periodic_table[atom.Z]['principal_quantum_number']**2)


def valence_state_indicator(atom):
    return vertex_degree(atom) + valence_degree(atom)


def kupchik_vertex_degree(atom):
    a = float(periodic_table[6]['covalent_radius'])/periodic_table[atom.Z]['covalent_radius']
    b = valence_electrones(atom) - sum([atom.Z == 1 for atom in atom.connected_with()])
    return a*b

def perturbation_delta_value(atom, perturbation_func):
    return valence_degree(atom) + sum([perturbation_func(atom, atom_n)*valence_degree(atom_n) for atom_n in atom.connected_with() if atom_n.Z != 1])

def intrinsic_state(atom):
    a = (2.0/periodic_table[atom.Z]["principal_quantum_number"]) **2
    b = a * valence_degree(atom) + 1
    c = vertex_degree(atom)
    return b/c


def madan_chemical_degree(atom):
    return sum([periodic_table[atom_n.Z]['relative_atomic_weight']/periodic_table[6]['relative_atomic_weight']
                for atom_n in atom.connected_with() if atom_n.Z != 1])


def extended_madan_vertex_degree(atom):
    s = 0
    for bond in atom.bonds:
        atom_n = filter(lambda x: x is not atom , [atom for atom in bond])[0]
        s += bond.conventional_bond_order*periodic_table[atom_n.Z]['relative_atomic_weight']/periodic_table[6]['relative_atomic_weight']
    return s

def chemical_extended_connectivity(atom):
    return sum([madan_chemical_degree(atom) for atom in atom.connected_with()])

def hu_xu_vertex_degree(atom):
    return vertex_degree(atom) * math.sqrt(atom.Z)

def consecutive_at_number(atom):
    return sum([math.sqrt(sqrt_2 + atom.Z) for atom in atom.connected_with()]) **2

def alikhanidi_vertex_degree(atom):
     return vertex_degree(atom) * consecutive_at_number(atom)

def augmented_valence(atom):
    """
    TODO Maybe not in this module received it from distance matrix it's characteristic of all molecule
    :param atom:
    :return:
    """

def ren_vertex_degree(atom):
    return vertex_degree(atom) + 1.0/(intrinsic_state(atom) * vertex_degree(atom))

def li_valence_vertex_degree(atom):
    return float(valence_electrones(atom) * valence_degree(atom))/(periodic_table[atom.Z]['principal_quantum_number']**2)

def yang_vertex_degree(atom):
    # atom in full molecule with hydrogen
    a = atom.Z - sum([atom.Z == 1 for atom in atom.connected_with()])
    b = vertex_degree(atom)
    c = periodic_table[atom.Z]['principal_quantum_number']**2 * periodic_table[atom.Z]['yang_electronegativity']
    return (float(a)*b)/c

def ct_vertex_degree(atom):
    a = vertex_degree(atom)
    b = reduce(operator.mul, [bond.conventional_bond_order for bond in atom.bonds])
    c = sum([(atom.Z + atom_n.Z)/12.0 for atom_n in atom.connected_with()])
    return math.sqrt(a*b*c)


def z_delta_number(atom):
    return float(valence_electrones(atom))/periodic_table[atom.Z]['principal_quantum_number']

def neignboor_interconenctivity(atom):
    descriptor =0
    neignboor = atom.connected_with()
    for atom1 in neignboor:
        descriptor += len(set(atom1.connected_with()) & set(neignboor))
    return descriptor/2

def cluster_coefficient_vertex(atom):
    v = vertex_degree(atom)
    if v > 1:
        return 2.0*neignboor_interconenctivity(atom)/(v* (v - 1))
    else:
        return 0


def edge_degree(bond):
    invariant = 0.0
    for atom in bond:
        invariant += vertex_degree(atom)
    return invariant - 2