# -*- coding: utf-8 -*-
"""
Can't classify category
"""
from collections import defaultdict

from utils.functional import cached
from vertex_degree import vertex_degree


@cached
def total_adjacency_index(molecule):
    return 2 * len(molecule.bonds)


@cached
def average_total_adjacency_index(molecule):
    return float(total_adjacency_index(molecule))/molecule.size


@cached
def density_index(molecule):
    return average_total_adjacency_index(molecule)/molecule.size


@cached
def no2_group_count(molecule):
    descriptor = 0
    for atom in molecule.atoms:
        if atom.Z == 7:
            connected_oxygens = filter(lambda atom: atom.Z == 8,  atom.connected_with())
            if len(connected_oxygens) == 2:
                descriptor += 1
    return descriptor


@cached
def etylenyl_group_count(molecule):
    descriptor = 0
    for bond in molecule.bonds:
        if bond.order == 2:
            atoms = [atom for atom in bond]
            if all (map(lambda x: x.Z ==6 , atoms)):
                feature = True
                for atom in atoms:
                    b = [bond_.order for bond_ in atom.bonds]
                    b.sort()
                    if b != [1,1,2]:
                        feature = False
                if feature:
                    descriptor += 1
    return descriptor


@cached
def hydroxi_group_count(molecule):
    descriptor = 0
    for atom in molecule.atoms:
        if atom.Z == 8:
            atoms = [atom.Z for atom in atom.connected_with()]
            if 1 in atoms:
                descriptor += 1
    return descriptor


@cached
def oxi_group_count(molecule):
    descriptor = 0
    for atom in molecule.atoms:
        if atom.Z == 8:
            if vertex_degree(atom) == 2:
                atoms = [atom.Z for atom in atom.connected_with()]
                if 1 not in atoms:
                    descriptor += 1
    return descriptor


@cached
def aldehyde_group_count(molecule):
    descriptor = 0
    for atom in molecule.atoms:
        h = False
        o = False
        if atom.Z == 6:
            for bond in atom.bonds:
                if bond.order == 2:
                    if 8 in [atom.Z for atom in bond]:
                        o = True
                if 1 in  [atom.Z for atom in bond]:
                    h = True
            if h and o:
                descriptor += 1

    return descriptor


@cached
def amino_groups_count(molecule):
    descriptor = 0
    for atom in molecule.atoms:
        if atom.Z == 7:
            if all(atom.Z in [1,6] for atom in atom.connected_with()):
                descriptor += 1
    return descriptor


@cached
def esters_group_count(molecule):
    descriptor = 0
    for atom in molecule.atoms:
        if atom.Z == 8:
            if vertex_degree(atom) == 2:
                atoms = atom.connected_with()
                if all([atom.Z == 6 for atom in atoms]):
                    for atom in atoms:
                        for bond in atom.bonds:
                            if bond.order ==2:
                                if any([atom.Z ==8 for atom in bond]):
                                    descriptor += 1
                                    continue
    return descriptor


@cached
def cetone_group_count(molecule):
    descriptor = 0
    for atom in molecule.atoms:
        if atom.Z == 8:
            if len(atom.bonds) == 1:
                c = atom.connected_with()[0]
                if c.Z ==6:
                    feature = True
                    for bond in c.bonds:
                        if bond.order == 1:
                            if not all([atom.Z != 1 and atom.Z != 8 for atom in bond]):
                                feature = False
                    if feature:
                        descriptor += 1
    return descriptor


@cached
def carboxy_group_count(molecule):
    descriptor = 0
    for atom in molecule.atoms:
        if atom.Z == 6:
            double_o = False
            hydroxi = False
            for bond in atom.bonds:
                if bond.order == 2:
                    if any([atom.Z == 8 for atom in bond]):
                        double_o = True
                if bond.order == 1:
                    if any([atom.Z == 8 for atom in bond]):
                        o_atom = filter(lambda x: x.Z == 8, [atom for atom in bond])[0]
                        if any([atom.Z == 1 for atom in o_atom.connected_with()]):
                            hydroxi = True
            if double_o and hydroxi:
                descriptor += 1
    return descriptor


@cached
def adsorbability_index(molecule):
    """
    Need article need system of priority

    :param molecule:
    :return:
    """
    descriptor = 0
    for atom in molecule.atoms:
        if atom.Z == 6:
            descriptor += 0.26
        elif atom.Z == 1:
            descriptor += 0.12
        elif atom.Z == 7:
            descriptor += 0.26
        elif atom.Z == 8:
            descriptor += 0.17
        elif atom.Z == 16:
            descriptor += 0.54
        elif atom.Z == 17:
            descriptor += 0.59
        elif atom.Z == 35:
            descriptor += 0.86
    descriptor += no2_group_count(molecule) * 0.21
    descriptor += etylenyl_group_count(molecule) * 0.19
    descriptor += hydroxi_group_count(molecule) * (-0.53)
    descriptor += oxi_group_count(molecule) * (-0.36)
    descriptor += aldehyde_group_count(molecule) * (-0.25)
    descriptor += amino_groups_count(molecule) * (-0.58)
    descriptor += esters_group_count(molecule) * (-0.28)
    descriptor += cetone_group_count(molecule) * (-0.30)
    descriptor += carboxy_group_count(molecule) * (-0.03)

    # unfinished
    return False