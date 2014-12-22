# -*- coding: utf-8 -*-
"""
Simple descriptor of molecule

"""
from utils.periodic_table import periodic_table
import math


def sv(molecule):
    return sum([(4 * math.pi * periodic_table[atom.Z]['vdw_radius']**3)/3 for atom in molecule.atoms])/((4 * math.pi * periodic_table[6]['vdw_radius']**3)/3)

def sp(molecule):
    return sum([periodic_table[atom.Z]['atomic_polarizability'] for atom in molecule.atoms])/periodic_table[6]['atomic_polarizability']

_ = {

    #molecular weight
    "MW": lambda x: x.molecular_mass(),

    # average molecular weight
    "AMW": lambda x: x.molecular_mass()/len(x.atoms),

    "SV": sv,

    "SP": sp
}