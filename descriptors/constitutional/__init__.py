# -*- coding: utf-8 -*-
"""
Constitutional descriptor (0D - 1D)
Molecular weight
Average molecular weight
Number of atoms
Number of bonds
Number of non hydrogen atoms
Number of non hydrogen bonds
Number of single bonds
Number of double bonds
Number of triple bonds
Number of aromatic bonds
Sum of conventional bond order
Rotatable bond count (J. Med. Chem. 2002, 45, 2615-2623)
Rotatable bond fraction
Number of H atoms
Number of B atoms
Number of C atoms
Number of N atoms
Number of O atoms
Number of P atoms
Number of S atoms
Number of Cl atoms
Number of Br atoms
Number of I atoms
Number of heavy atoms
Number of heteroatoms
Number of halogen atoms
Percentage of H atoms
Percentage of C atoms
Percentage of N atoms
Percentage of O atoms
Percentage of halogen atoms
Number of C-sp3 atoms
Number of C-sp2 atoms
Number of C-sp atoms
Sum of van der Waals Volumes (scaled on Carbon atom)
Mean of van der Waals Volumes (scaled on Carbon atom)
Sum of atomic polarizability (scaled on Carbon atom)
Mean of atomic polarizability (scaled on Carbon atom)
Sum of first_ionization_potential (scaled on Carbon atom)
Mean of first_ionization_potential (scaled on Carbon atom)
Sum of Sanderson electronegativity (scaled on Carbon atom)
Mean of Sanderson electronegativity (scaled on Carbon atom)
"""
from utils.periodic_table import periodic_table
import math
from calc.graph import cycle_bonds
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

def rotatable_bond_count(molecule):
    """
    number of rotatable bonds
    :param molecule:
    :return:
    """
    descriptor = 0
    molecule = molecule.hydrogen_suppressed
    cycle_bond_list = cycle_bonds(molecule)
    for bond in molecule.bonds:
        rotatable = True
        if bond.order > 1: # double, triple bond not count
            rotatable = False
        elif bond in cycle_bond_list:
            rotatable = False
        elif bond.is_amide_bond():
            rotatable = False
        elif bond.is_terminal():
            rotatable = False
        elif bond.adjacent_to_triple_bond():
            rotatable = False
        if rotatable:
            descriptor += 1
    return descriptor

def rotatable_bond_fraction(molecule):
    """
    rotatable bond fraction
    :param molecule:
    :return:
    """
    return float(rotatable_bond_count(molecule))/molecule.size

def number_of_hydrogen_atoms(molecule):
    return len(filter( lambda x: x.Z == 1 ,molecule.atoms))

def number_of_boron_atoms(molecule):
    return len(filter( lambda x: x.Z == 5 ,molecule.atoms))

def number_of_carbon_atoms(molecule):
    return len(filter( lambda x: x.Z == 6 ,molecule.atoms))

def number_of_nytrogen_atoms(molecule):
    return len(filter( lambda x: x.Z == 7,molecule.atoms))

def number_of_oxygen_atoms(molecule):
    return len(filter( lambda x: x.Z == 8 ,molecule.atoms))

def number_of_phosphorous_atoms(molecule):
    return len(filter( lambda x: x.Z == 15 ,molecule.atoms))

def number_of_sulfur_atoms(molecule):
    return len(filter( lambda x: x.Z == 16 ,molecule.atoms))

def number_of_chlorine_atoms(molecule):
    return len(filter( lambda x: x.Z == 17 ,molecule.atoms))

def number_of_bromine_atoms(molecule):
    return len(filter( lambda x: x.Z == 35 ,molecule.atoms))

def number_of_iodine_atoms(molecule):
    return len(filter( lambda x: x.Z == 53 ,molecule.atoms))

def number_of_heavy_atoms(molecule):
    return len(filter( lambda x: x.Z != 1 ,molecule.atoms))

def number_of_heteroatoms(molecule):
    return len(filter( lambda x: x.Z  not in [1, 6] ,molecule.atoms))

def number_of_halogen_atoms(molecule):
    return len(filter( lambda x: x.Z in [9, 17, 35, 53] ,molecule.atoms))

def percentage_of_hydrogen_atoms(molecule):
    return 100*float(number_of_hydrogen_atoms(molecule))/molecule.size

def percentage_of_carbon_atoms(molecule):
    return 100*float(number_of_carbon_atoms(molecule))/molecule.size

def percentage_of_nytrogen_atoms(molecule):
    return 100*float(number_of_nytrogen_atoms(molecule))/molecule.size

def percentage_of_oxygen_atoms(molecule):
    return 100*float(number_of_oxygen_atoms(molecule))/molecule.size

def percentage_of_halogen(molecule):
    return 100*float(number_of_halogen_atoms(molecule))/molecule.size

def number_of_csp3_carbon_atoms(molecule):
    descriptor = 0
    for atom in molecule.atoms:
        if atom.Z != 6 or atom.aromatic:
            continue
        not_sp3 = False
        for bond in atom.bonds:
            if bond.order != 1:
                not_sp3 = True
                break
        descriptor += 1 if not bad else 0
    return descriptor

def number_of_csp2_carbon_atoms(molecule):
    descriptor = 0
    for atom in molecule.atoms:
        if atom.Z != 6:
            continue
        if atom.aromatic:
            descriptor += 1
            continue
        if any([bond.order == 2 for bond in atom.bonds]):
            descriptor += 1
    return descriptor

def number_of_csp_carbon_atoms(molecule):
    descriptor = 0
    for atom in molecule.atoms:
        if atom.Z != 6 or atom.aromatic:
            continue
        if any([bond.order == 3 for bond in atom.bonds]):
            descriptor += 1
    return descriptor
