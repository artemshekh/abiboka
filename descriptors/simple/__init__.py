# -*- coding: utf-8 -*-
"""
Simple descriptor of molecule

"""
from utils.periodic_table import periodic_table
import math


def sv(molecule):
    return sum([(4 * math.pi * periodic_table[atom.Z]['vdw_radius']**3)/3 for atom in molecule.atoms])/((4 * math.pi * periodic_table[6]['vdw_radius']**3)/3)

def mv(molecule):
    return sv(molecule)/len(molecule.atoms)

def sp(molecule):
    return sum([periodic_table[atom.Z]['atomic_polarizability'] for atom in molecule.atoms])/periodic_table[6]['atomic_polarizability']

def mp(molecule):
    return sp(molecule)/len(molecule.atoms)

def si(molecule):
    return sum([periodic_table[atom.Z]['1st_ionization_energy'] for atom in molecule.atoms])/periodic_table[6]['1st_ionization_energy']

def mi(molecule):
    return si(molecule)/len(molecule.atoms)

def nbo(molecule):
    nbo_ = 0
    for bond in molecule.bonds:
        h_in_bond = False
        for atom in bond:
            if atom.Z == 1:
                h_in_bond = True
        if not h_in_bond:
            nbo_ += 1
    return nbo_

def scbo(molecule):
    scbo_ = 0
    for bond in molecule.bonds:
        if bond.is_aromatic():
            scbo_ += 1.5
        else:
            scbo_ += bond.order
    return scbo_

_ = {

    #molecular weight
    "MW": lambda x: x.molecular_mass(),

    # average molecular weight
    "AMW": lambda x: x.molecular_mass()/len(x.atoms),

    "SV": sv,

    "MV": mv,

    "SP": sp,

    "MP": mp,

    "SI": si,

    "MI": mi,

    # number of atoms
    "nAT": lambda x: len(x.atoms),

    "nSK": lambda x: len(filter( lambda x: x.Z !=1 ,x.atoms)),

    "nBT": lambda x: len(x.bonds),

    "nBO": nbo,

    # number of multiple bonds (2 or 3) aromatic bonds is not multiple. Good question?
    "nBM": lambda x: len(filter( lambda x: x.order !=1 ,x.bonds)),

    "SCBO": scbo,

}