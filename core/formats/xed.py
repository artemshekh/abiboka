# -*- coding: utf-8 -*-
"""
http://www.cresset-group.com/products/xedtools/
File format definition: http://www.cresset-bmd.com/fieldscreen_manual/File_format.pdf
Atom types: http://www.cresset-bmd.com/fieldscreen_manual/Xed_atom_types.pdf
# TODO aromacity n,o,s !!!!!!!!!
"""

import re

from structure.Molecule import Molecule
from structure.Atom import Atom
from structure.Bond import Bond
from parser import Parser


class XedParser(Parser):

    def __init__(self):
        self.CONNECTIVITY_STRING = re.compile('^(\ *\d*\ *\d)*$')
        self.bond_order_converter = {
            (1, 1): 1,
(1, 2): 2,
(1, 3): 1,
(1, 4): 1,
(1, 5): 1,
(1, 6): 1,
(1, 8): 1,
(1, 9): 1,
(1, 10): 1,
(1, 12): 1,
(1, 14): 1,
(1, 15): 1,
(1, 16): 1,
(1, 17): 1,
(1, 18): 1,
(1, 19): 1,
(1, 23): 1,
(1, 25): 1,
(1, 26): 1,
(1, 37): 1,
(1, 40): 1,
(2, 2): 2,
(2, 3): 1,
(2, 4): 1,
(2, 7): 2,
(2, 9): 1,
(2, 10): 1,
(2, 11): 2,
(2, 12): 1,
(2, 13): 2,
(2, 15): 1,
(2, 16): 1,
(2, 17): 1,
(2, 18): 1,
(2, 19): 1,
(2, 25): 1,
(2, 26): 1,
(3, 3): 1,
(3, 4): 1,
(3, 5): 1,
(3, 8): 1,
(3, 9): 1,
(3, 10): 1,
(3, 12): 1,
(3, 14): 2,
(3, 15): 1,
(3, 16): 1,
(3, 17): 1,
(3, 18): 1,
(3, 19): 1,
(3, 23): 1,
(3, 25): 1,
(3, 26): 2,
(3, 40): 1,
(4, 4): 3,
(4, 8): 2,
(4, 9): 1,
(4, 10): 3,
(4, 11): 2,
(4, 12): 1,
(4, 13): 2,
(4, 14): 1,
(4, 15): 1,
(4, 18): 1,
(4, 23): 3,
(4, 25): 1,
(4, 26): 1,
(5, 12): 1,
(5, 16): 1,
(5, 26): 1,
(6, 6): 1,
(6, 14): 1,
(6, 15): 1,
(7, 23): 2,
(8, 9): 1,
(8, 10): 1,
(8, 11): 2,
(8, 12): 1,
(8, 15): 1,
(8, 17): 1,
(8, 25): 1,
(9, 9): 1,
(9, 10): 1,
(9, 12): 1,
(9, 14): 1,
(9, 15): 1,
(9, 16): 1,
(9, 17): 1,
(9, 18): 1,
(9, 19): 1,
(9, 22): 2,
(9, 25): 1,
(9, 26): 1,
(10, 10): 1,
(10, 12): 1,
(10, 14): 1,
(10, 15): 1,
(10, 25): 1,
(10, 37): 1,
(10, 40): 1,
(12, 12): 1,
(12, 14): 1,
(12, 15): 1,
(12, 16): 1,
(12, 17): 1,
(12, 18): 1,
(12, 24): 2,
(13, 14): 2,
(14, 14): 1,
(14, 15): 1,
(14, 17): 1,
(14, 18): 1,
(14, 19): 1,
(14, 23): 2,
(14, 24): 2,
(16, 25): 1,
(23, 23): 2,
(25, 26): 1,
        }

    def decode(self, file):
        connectivity_list = []  # list of bonds in molecule [(1,3), (2, 7)]
        molecule = Molecule()
        atom_count = 1
        atom_dct = {}
        atom_types = {}

        for line in file:
            if re.match(self.CONNECTIVITY_STRING, line):  # string containing connectivity information
                bonds = map(lambda x: int(x), filter(lambda x: x != '', line.strip('\n').split(' ')))
                bonds = zip(*[iter(bonds)]*2)
                connectivity_list += bonds

            if 'MOL' in line and 'X' not in line:
                atom_line = filter(lambda x: x != '', line.strip('\n').split(' '))
                z = atom_line[0]
                x_coord, y_coord, z_coord = atom_line[1], atom_line[2], atom_line[3]
                atom_type = atom_line[4]

                aromatic = False
                if atom_type == '3':
                    aromatic = True

                atom = Atom(z=int(z), aromatic=aromatic)
                atom.coords = (float(x_coord), float(y_coord), float(z_coord))
                molecule.add_atom(atom)

                atom_dct[atom_count] = atom
                atom_types[atom_count] = int(atom_type)
                atom_count += 1

        size = molecule.size
        for bond in connectivity_list:
            if bond[0] < size and bond[1] < size:
                atom1 = atom_dct[bond[0]]
                atom2 = atom_dct[bond[1]]
                key = [atom_types[bond[0]], atom_types[bond[1]]]
                key.sort()
                key = tuple(key)
                if key in self.bond_order_converter:
                    order = self.bond_order_converter[key]
                else:
                    order = 0  # ubdefined order

                b = Bond(order=order, atom1=atom1, atom2=atom2)
                molecule.bonds.append(b)
        for bond in molecule.bonds:
            atoms = [atom.Z for atom in bond]
            if not bond.order:
                if 6 in atoms:
                    for atom in bond:
                        if atom.Z == 6:
                            v = 0
                            for b in atom.bonds:
                                v += b.order
                            bond.order = 4 - v or 1
                            break
                elif 7 in atoms:
                    for atom in bond:
                        if atom.Z == 7:
                            v = 0
                            for b in atom.bonds:
                                v += b.order
                            if v < 3:
                                bond.order = 3-v
                            else:
                                bond.order = 5-v
                            break
        return molecule