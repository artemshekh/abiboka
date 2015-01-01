# -*- coding: utf-8 -*-
"""
Simple descriptor of molecule

"""
from utils.periodic_table import periodic_table
import math
from calc.graph import check_cycle, dfs


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
    for bond in molecule.hydrogen_suppressed.bonds:
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
    molecule = molecule.hydrogen_suppressed
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

def nab(molecule):
    _ = 0
    for bond in molecule.bonds:
        _ += 1 if bond.is_aromatic() else 0
    return _

def ncsp3(molecule):
    _ = 0
    for atom in molecule.atoms:
        if atom.Z != 6 or atom.aromatic:
            continue
        bad = False
        for bond in atom.bonds:
            if bond.order != 1:
                bad = True
                break
        _ += 1 if not bad else 0
    return _

def ncsp2(molecule):
    _ = 0
    for atom in molecule.atoms:
        if atom.Z != 6:
            continue
        if atom.aromatic:
            _ += 1
            continue
        if any([bond.order == 2 for bond in atom.bonds]):
            _ += 1
    return _

def ncsp(molecule):
    _ = 0
    for atom in molecule.atoms:
        if atom.Z != 6 or atom.aromatic:
            continue
        if any([bond.order == 3 for bond in atom.bonds]):
            _ += 1
    return _




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

    # number of non-H atoms
    "nSK": lambda x: len(filter( lambda x: x.Z !=1 ,x.atoms)),

    "nBT": lambda x: len(x.bonds),

    "nBO": nbo,

    # number of multiple bonds (2 or 3) aromatic bonds is not multiple. Good question?
    "nBM": lambda x: len(filter( lambda x: x.order !=1 ,x.bonds)),

    "SCBO": scbo,

    "RBN": rbn,

    "RBF": rbf,

    # number of double bonds
    "NDB": lambda x: len([bond.order for bond in x.bonds if bond.order==2]),

    # number of triple bonds
    "NTB": lambda x: len([bond.order for bond in x.bonds if bond.order==3]),

    "NAB": nab,

    # number of Hydrogen atoms
    "nH": lambda x: len(filter( lambda x: x.Z ==1 ,x.atoms)),

    # number of Carbon atoms
    "nC": lambda x: len(filter( lambda x: x.Z ==6 ,x.atoms)),

    # number of Nytrogen atoms
    "nN": lambda x: len(filter( lambda x: x.Z ==7 ,x.atoms)),

    # number of Oxygen atoms
    "nO": lambda x: len(filter( lambda x: x.Z ==8 ,x.atoms)),

    # number of Posphorus atoms
    "nP": lambda x: len(filter( lambda x: x.Z ==15 ,x.atoms)),

    # number of Sulfur atoms
    "NS": lambda x: len(filter( lambda x: x.Z ==16 ,x.atoms)),

    # number of Fluorine atoms
    "nF": lambda x: len(filter( lambda x: x.Z ==9 ,x.atoms)),

    # number of Chlorine atoms
    "nCL": lambda x: len(filter( lambda x: x.Z ==17 ,x.atoms)),

    # number of Bromine atoms
    "nBR": lambda x: len(filter( lambda x: x.Z ==35 ,x.atoms)),

    # number of Iodine atoms
    "nI": lambda x: len(filter( lambda x: x.Z ==53 ,x.atoms)),

    # number of Boron atoms
    "nB": lambda x: len(filter( lambda x: x.Z ==5 ,x.atoms)),

    # number of heavy atom ( see nSK)
    "nHM": lambda x: len(filter( lambda x: x.Z !=1 ,x.atoms)),

    # number of heteroatoms
    "nHet": lambda x: len(filter(lambda x: x.Z not in [1,6], x.atoms)),

    # number of halogen atoms
    "nX": lambda x: len(filter(lambda x: x.Z  in [9,17,35,53], x.atoms)),

    # percentage of H atoms
    "H%": lambda x: 100*float(len(filter( lambda x: x.Z ==1 ,x.atoms)))/len(x.atoms),

    # percentage of C atoms
    "C%": lambda x: 100*float(len(filter( lambda x: x.Z ==6 ,x.atoms)))/len(x.atoms),

    # percentage of N atoms
    "N%": lambda x: 100*float(len(filter( lambda x: x.Z ==7 ,x.atoms)))/len(x.atoms),

    # percentage of O atoms
    "O%": lambda x: 100*float(len(filter( lambda x: x.Z ==8 ,x.atoms)))/len(x.atoms),

    # percentage of halogen atoms
    "X%": lambda x: 100*float(len(filter( lambda x: x.Z in [9,17,35,53] ,x.atoms)))/len(x.atoms),

    # number of Csp3 Carbon atoms
    "nCsp3": ncsp3,

    # number of Csp2 Carbon atoms
    "nCsp2": ncsp2,

    # number of Csp Carbon atoms
    "nCsp": ncsp,

    # cyclomatic number
    "nCIC": cyclomatic_number,

}