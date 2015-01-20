# -*- coding: utf-8 -*-
"""
Molecule for testing purposes
"""
from structure.Atom import Atom
from structure.Bond import Bond
from structure.Molecule import Molecule


def ethanol():
    m = Molecule()
    m.atoms += [Atom(z=6) for x in range(2)]
    m.atoms += [Atom(z=1) for x in range(6)]
    m.atoms += [Atom(z=8) for x in range(1)]
    atoms = m.atoms
    bonds = [(0, 1), (0, 2), (0, 3), (0, 4), (1, 8), (1, 5), (1, 6), (7, 8)]
    for i, j in bonds:
        bond = Bond(atom1=atoms[i], atom2=atoms[j])
        m.bonds.append(bond)
    return m


def ethynylbenzoicacid():
    m = Molecule()
    m.atoms += [Atom(z=6) for x in range(9)]
    m.atoms += [Atom(z=1) for x in range(6)]
    m.atoms += [Atom(z=8) for x in range(2)]
    atoms = m.atoms
    aromatic_bonds = [(2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (2, 7)]
    for i, j in aromatic_bonds:
        atoms[i].aromatic = True
        atoms[j].aromatic = True
        bond = Bond(atom1=atoms[i], atom2=atoms[j])
        m.bonds.append(bond)
    single_bonds = [(0, 9), (1, 2), (3, 10), (4, 11), (5, 12), (7, 13), (6, 8), (8, 16), (14, 16)]
    for i, j in single_bonds:
        bond = Bond(atom1=atoms[i], atom2=atoms[j])
        m.bonds.append(bond)
    double_bond = [(8, 15)]
    for i, j in double_bond:
        bond = Bond(order=2, atom1=atoms[i], atom2=atoms[j])
        m.bonds.append(bond)
    triple_bond = [(0, 1)]
    for i, j in triple_bond:
        bond = Bond(order=3, atom1=atoms[i], atom2=atoms[j])
        m.bonds.append(bond)
    return m


def borane():
    m = Molecule()
    m.atoms += [Atom(z=5) for x in range(1)]
    m.atoms += [Atom(z=1) for x in range(3)]
    bonds = [(0, 1), (0, 2), (0, 3)]
    atoms = m.atoms
    for i, j in bonds:
        bond = Bond(atom1=atoms[i], atom2=atoms[j])
        m.bonds.append(bond)
    return m


def ammiac():
    m = Molecule()
    m.atoms += [Atom(z=7) for x in range(1)]
    m.atoms += [Atom(z=1) for x in range(3)]
    bonds = [(0, 1), (0, 2), (0, 3)]
    atoms = m.atoms
    for i, j in bonds:
        bond = Bond(atom1=atoms[i], atom2=atoms[j])
        m.bonds.append(bond)
    return m


def phosphane():
    m = Molecule()
    m.atoms += [Atom(z=15) for x in range(1)]
    m.atoms += [Atom(z=1) for x in range(3)]
    bonds = [(0, 1), (0, 2), (0, 3)]
    atoms = m.atoms
    for i, j in bonds:
        bond = Bond(atom1=atoms[i], atom2=atoms[j])
        m.bonds.append(bond)
    return m


def hydrogen_sulfide():
    m = Molecule()
    m.atoms += [Atom(z=16) for x in range(1)]
    m.atoms += [Atom(z=1) for x in range(2)]
    bonds = [(0, 1), (0, 2)]
    atoms = m.atoms
    for i, j in bonds:
        bond = Bond(atom1=atoms[i], atom2=atoms[j])
        m.bonds.append(bond)
    return m


def hydrogen_chloride():
    m = Molecule()
    m.atoms += [Atom(z=17) for x in range(1)]
    m.atoms += [Atom(z=1) for x in range(1)]
    bonds = [(0, 1)]
    atoms = m.atoms
    for i, j in bonds:
        bond = Bond(atom1=atoms[i], atom2=atoms[j])
        m.bonds.append(bond)
    return m


def hydrogen_bromide():
    m = Molecule()
    m.atoms += [Atom(z=35) for x in range(1)]
    m.atoms += [Atom(z=1) for x in range(1)]
    bonds = [(0, 1)]
    atoms = m.atoms
    for i, j in bonds:
        bond = Bond(atom1=atoms[i], atom2=atoms[j])
        m.bonds.append(bond)
    return m


def hydrogen_iodide():
    m = Molecule()
    m.atoms += [Atom(z=53) for x in range(1)]
    m.atoms += [Atom(z=1) for x in range(1)]
    bonds = [(0, 1)]
    atoms = m.atoms
    for i, j in bonds:
        bond = Bond(atom1=atoms[i], atom2=atoms[j])
        m.bonds.append(bond)
    return m

