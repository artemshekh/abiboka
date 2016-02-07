# -*- coding: utf-8 -*-
"""
Schrodinger mae format

"""
import re

from parser import Parser
from structure.Molecule import Molecule
from structure.Atom import Atom
from structure.Bond import Bond


class MaeParser(Parser):

    def __init__(self):
        self.OPEN_CURLY_BRACKET = re.compile('.*\{[ ]*\\n')
        self.CLOSED_CURLY_BRACKET = '} \n'

    def encode(self, molecule):
        raise NotImplementedError

    def decode(self, filename):
        mae_file = open(filename, 'r')
        current_block = Block('root', None)

        for line in mae_file:
            if line == '\n':
                continue
            if self.OPEN_CURLY_BRACKET.match(line):
                block = self.create_block(line, current_block)
                current_block.blocks.append(block)
                current_block = block
            elif self.CLOSED_CURLY_BRACKET in line:
                # parse block
                # return molecule if its f_m_ct_block
                current_block = current_block.parent_block
                if current_block.name == 'root':
                    last_block = current_block.blocks[-1]
                    if last_block.name == 'f_m_ct':
                        mol = self.parse_fmct_block(last_block)
                        yield mol
            else:
                block.content.append(line)

    def create_block(self, line, parent_block):
        name = line.split('{')[0].strip(' ')
        block = Block(name, parent_block)
        return block

    def parse_fmct_block(self, fmct_block):
        mol = Molecule()
        atoms_by_index = {}
        for block in fmct_block.blocks:
            if block.name.startswith('m_atom'):
                columns_by_index = {'atom_index': 0}
                parse_names = True
                count = 1
                for line in block.content:
                    line = line.strip('\n').strip(' ')
                    if line.startswith('#'):
                        continue
                    if ':::' == line:
                        parse_names = False
                        continue
                    if parse_names:
                        columns_by_index[line] = count
                        count += 1
                    else:
                        line = line.replace('"', '')
                        data = [x for x in line.split(' ') if x != '']
                        atomic_number = int(data[columns_by_index['i_m_atomic_number']])
                        x_coord = float(data[columns_by_index['r_m_x_coord']])
                        y_coord = float(data[columns_by_index['r_m_y_coord']])
                        z_coord = float(data[columns_by_index['r_m_z_coord']])
                        atom = Atom(z=atomic_number)
                        atom.x, atom.y, atom.z = x_coord, y_coord, z_coord
                        atom.mae_props = {}
                        for column, index in columns_by_index.items():
                            atom.mae_props[column] = data[index]
                        mol.atoms.append(atom)
                        atoms_by_index[int(data[0])] = atom
            elif block.name.startswith('m_bond'):
                columns_by_index = {'bond_index': 0}
                parse_names = True
                count = 1
                for line in block.content:
                    line = line.strip('\n').strip(' ')
                    if line.startswith('#'):
                        continue
                    if ':::' == line:
                        parse_names = False
                        continue
                    if parse_names:
                        columns_by_index[line] = count
                        count += 1
                    else:
                        if not line:
                            continue
                        line = line.replace('"', '')
                        data = [x for x in line.split(' ') if x != '']
                        atom1 = atoms_by_index[int(data[int(columns_by_index['i_m_from'])])]
                        atom2 = atoms_by_index[int(data[int(columns_by_index['i_m_to'])])]
                        bond_order = int(data[int(columns_by_index['i_m_order'])])
                        b = Bond(order=bond_order, atom1=atom1, atom2=atom2)
                        mol.bonds.append(b)
        return mol


class Block(object):

    def __init__(self, name, parent_block):
        self.name = name
        self.parent_block = parent_block
        self.content = []
        self.blocks = []