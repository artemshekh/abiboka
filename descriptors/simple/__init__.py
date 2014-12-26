# -*- coding: utf-8 -*-
"""
Simple descriptor of molecule

"""
from utils.periodic_table import periodic_table
import math
from calc.graph import check_cycle


def sv(molecule):
    """
    Sum of Van der Waals volume
    :param molecule:
    :return:
    """
    return sum([(4 * math.pi * periodic_table[atom.Z]['vdw_radius']**3)/3 for atom in molecule.atoms])/((4 * math.pi * periodic_table[6]['vdw_radius']**3)/3)

def mv(molecule):
    """
    Mean of Van der Waals volume
    :param molecule:
    :return:
    """
    return sv(molecule)/len(molecule.atoms)

def sp(molecule):
    """
    sum of atomic polarizabilities (scaled on Carbon atom)
    :param molecule:
    :return:
    """
    return sum([periodic_table[atom.Z]['atomic_polarizability'] for atom in molecule.atoms])/periodic_table[6]['atomic_polarizability']

def mp(molecule):
    """
    mean atomic polarizability (scaled on Carbon atom)
    :param molecule:
    :return:
    """
    return sp(molecule)/len(molecule.atoms)

def si(molecule):
    """
    sum of first ionization potentials (scaled on Carbon atom)
    :param molecule:
    :return:
    """
    return sum([periodic_table[atom.Z]['1st_ionization_energy'] for atom in molecule.atoms])/periodic_table[6]['1st_ionization_energy']

def mi(molecule):
    """
    mean first ionization potential (scaled on Carbon atom)
    :param molecule:
    :return:
    """
    return si(molecule)/len(molecule.atoms)

def nbo(molecule):
    """
    number of non-H bonds
    :param molecule:
    :return:
    """
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
    """
    sum of conventional bond orders (H-depleted)
    :param molecule:
    :return:
    """
    scbo_ = 0
    for bond in molecule.hydrogen_suppressed().bonds:
        if bond.is_aromatic():
            scbo_ += 1.5
        else:
            scbo_ += bond.order
    return scbo_

def rbn(molecule):
    """
    number of rotatable bonds
    :param molecule:
    :return:
    """
    rbn_ = 0
    molecule = molecule.hydrogen_suppressed()
    cycle_atoms = check_cycle(molecule)
    for bond in molecule.bonds:
        rotatable = True
        if bond.order > 1:
            rotatable = False
        atoms = [atom for atom in bond]

        if atoms[0] in cycle_atoms and atoms[1] in cycle_atoms:
            rotatable = False
        elif atoms[0].Z == 7 or atoms[1].Z == 7:
            rotatable = False
        elif len(atoms[0].bonds) == 1 or len(atoms[1].bonds) == 1:
            rotatable = False
        if rotatable:
            rbn_ += 1
    return rbn_

def rbf(molecule):
    """
    rotatable bond fraction
    :param molecule:
    :return:
    """
    return float(rbn(molecule))/len(molecule.bonds)

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

    "RBN": rbn,

    "RBF": rbf,

}