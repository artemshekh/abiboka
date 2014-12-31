# -*- coding: utf-8 -*-
"""
Topological descriptors
"""
import operator
from utils.periodic_table import periodic_table


def zm1(molecule):
    return sum([len(atom.bonds)**2 for atom in molecule.hydrogen_suppressed().atoms])

def zm1_H(molecule):
    return sum([len(atom.bonds)**2 for atom in moleculeatoms])

close_shell = [1, 3, 11, 19, 37, 55, 87]

def valence_electrones(atom):

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


def zm1v(molecule):
    return sum([valence_degree(atom)**2 for atom in molecule.hydrogen_suppressed().atoms])

def zm1v_(molecule):
    # depends from core electrons
    return sum([valence_degree(atom)**2 for atom in molecule.hydrogen_suppressed().atoms])

def kupchik_degree(atom):
    return  (float(periodic_table[6]['covalent_radius'])/periodic_table[atom.Z]["covalent_radius"])* valence_degree(atom)

def zm1kup(molecule):
    return sum([kupchik_degree(atom)**2 for atom in molecule.atoms])

def madan_degree(atom):
    return sum([periodic_table[atom1.Z]["relative_atomic_weight"]/periodic_table[6]["relative_atomic_weight"] for atom1 in atom.connected_with()])

def zm1mad(molecule):
    return sum([madan_degree(atom) ** 2 for atom in molecule.atoms])

def permutation_additive(atom, perm_coefficient):
    return valence_degree(atom) + sum([perm_coefficient * valence_degree(atom) for _ in atom.connected_with()])

def permutation_multiplicative(atom, perm_coefficient):
    return valence_degree(atom) + reduce(operator.mul, [perm_coefficient * valence_degree(atom) for _ in atom.connected_with()])

def zm1per(molecule, permutation_coefficient):
    return sum([permutation_additive(atom, permutation_coefficient) ** 2 for atom in molecule.atoms])

def zm1mulper(molecule, permutation_coefficient):
    return sum([permutation_multiplicative(atom, permutation_coefficient) ** 2 for atom in molecule.atoms])

def zm2(molecule):
    _ = []
    for bond in molecule.hydrogen_suppressed().bonds:
        atoms = [atom for atom in bond]
        _.append(len(atoms[0].bonds)*len(atoms[1].bonds))
    return sum(_)

def zm2v(molecule):
    _ = []
    for bond in molecule.bonds:
        atoms = [atom for atom in bond]
        _.append(valence_degree(atoms[0])*valence_degree(atoms[1]))
    return sum(_)

def zm2kup(molecule):
    _ = []
    for bond in molecule.bonds:
        atoms = [atom for atom in bond]
        _.append(kupchik_degree(atoms[0])*kupchik_degree(atoms[1]))
    return sum(_)


def zm2mad(molecule):
    _ = []
    for bond in molecule.bonds:
        atoms = [atom for atom in bond]
        _.append(madan_degree(atoms[0])*madan_degree(atoms[1]))
    return sum(_)

def zm2per(molecule, permutative_coefficent):
    _ = []
    for bond in molecule.bonds:
        atoms = [atom for atom in bond]
        _.append(permutation_additive(atoms[0], permutative_coefficent)*permutation_additive(atoms[1], permutative_coefficent))
    return sum(_)

def zm2mulper(molecule, permutative_coefficent):
    _ = []
    for bond in molecule.bonds:
        atoms = [atom for atom in bond]
        _.append(permutation_multiplicative(atoms[0], permutative_coefficent)*permutation_multiplicative(atoms[1], permutative_coefficent))
    return sum(_)





















