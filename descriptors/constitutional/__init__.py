# -*- coding: utf-8 -*-
"""
Simple descriptor of molecule

"""
from utils.periodic_table import periodic_table
import math
from calc.graph import check_cycle
from descriptors.ring_descriptor import cyclomatic_number


def molecular_weight(molecule):
    return molecule.molecular_mass()

def average_molecular_weight(molecule):
    return molecular_weight(molecule)/molecule.size

def sum_van_der_waals_volume(molecule):
    """
    Sum of Van der Waals volume (scaled on Carbon atom)
    :param molecule:
    :return:
    """
    constant = (4.0/3)*math.pi
    volume_carbon_atom = constant*(periodic_table[6]['vdw_radius']**3)
    sum_volume = sum([constant * periodic_table[atom.Z]['vdw_radius']**3  for atom in molecule.atoms])
    return sum_volume/volume_carbon_atom

def mean_van_der_waals_volume(molecule):
    """
    Mean of Van der Waals volume
    :param molecule:
    :return:
    """
    return sum_van_der_waals_volume(molecule)/molecule.size

def sum_of_atom_polarizability(molecule):
    """
    sum of atomic polarizabilities (scaled on Carbon atom)
    :param molecule:
    :return:
    """
    polarizability_carbon_atom = periodic_table[6]['atomic_polarizability']
    sum_polarizability = sum([periodic_table[atom.Z]['atomic_polarizability'] for atom in molecule.atoms])
    return sum_polarizability/polarizability_carbon_atom

def mean_sum_atom_polarizability(molecule):
    """
    mean atomic polarizability (scaled on Carbon atom)
    :param molecule:
    :return:
    """
    return sum_of_atom_polarizability(molecule)/molecule.size

def sum_of_first_ionization_potentials(molecule):
    """
    sum of first ionization potentials (scaled on Carbon atom)
    :param molecule:
    :return:
    """
    carbon_first_ionization_potential = periodic_table[6]['1st_ionization_energy']
    sum_first_ionization_potential_molecule = sum([periodic_table[atom.Z]['1st_ionization_energy'] for atom in molecule.atoms])
    return sum_first_ionization_potential_molecule/carbon_first_ionization_potential

def mean_sum_of_first_ionization_potentials(molecule):
    """
    mean first ionization potential (scaled on Carbon atom)
    :param molecule:
    :return:
    """
    return sum_of_first_ionization_potentials(molecule)/molecule.size

def sum_of_sanderson_electronegativity(molecule):
    sanderson_electronegativity_c = periodic_table[6]['sanderson_electronegativity']
    sanderson_electronegativity_molecule = sum([periodic_table[atom.Z]['sanderson_electronegativity'] for atom in molecule.atoms])
    return sanderson_electronegativity_molecule/sanderson_electronegativity_c

def mean_sum_of_sanderson_electronegativity(molecule):
    return sum_of_sanderson_electronegativity(molecule)/molecule.size

def number_of_atoms(molecule):
    return molecule.size

def number_of_non_hydrogen_atoms(molecule):
    return len(filter( lambda atom: atom.Z != 1 , molecule.atoms))

def number_of_bonds(molecule):
    return len(molecule.bonds)

def number_of_non_hydrogen_bonds(molecule):
    return len(molecule.hydrogen_suppressed.bonds)

def number_of_single_bond(molecule):
    descriptor = 0
    for bond in molecule.bonds:
        if bond.order == 1 and not bond.is_aromatic():
            print bond
            descriptor += 1
    return descriptor

def number_of_double_bonds(molecule):
    descriptor = 0
    for bond in molecule.bonds:
        if bond.order == 2:
            descriptor += 1
    return descriptor

def number_of_triple_bonds(molecule):
    descriptor = 0
    for bond in molecule.bonds:
        if bond.order == 3:
            descriptor += 1
    return descriptor

def number_of_aromatic_bonds(molecule):
    descriptor = 0
    for bond in molecule.bonds:
        if bond.is_aromatic():
            descriptor += 1
    return descriptor

def sum_of_conventional_bond_order(molecule):
    descriptor = 0
    for bond in molecule.hydrogen_suppressed.bonds:
        if bond.is_aromatic():
            descriptor += 1.5
        else:
            descriptor += bond.order
    return descriptor

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