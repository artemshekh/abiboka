# -*- coding: utf-8 -*-
"""
Schrodinger mae format

"""
import os

from parser import Parser
from structure.Molecule import Molecule
from structure.Atom import Atom
from structure.Bond import Bond


class MaeParser(Parser):

    def __init__(self):
        self.header = open('templates/mae_header', 'r').read()


    def decode(self, file):
        block = ''
        level = 0
        for line in open(file, 'r'):
            #print line, level
            if 'f_m_ct' in line:
                # start of block
                level += 1
                block += line
            elif '{' in line and level != 0:
                level += 1
                block += line
            elif '}' in line and level != 0:
                level -= 1
                block += line
                if level == 0:
                    atoms_block, bonds_block = self.parse_f_m_ct(block)
                    atoms = self.parse_atoms_block(atoms_block)
                    bonds = self.parse_bond_block(bonds_block)
                    mol = self.construct_mol(atoms, bonds)
                    yield mol

                    block = ''
            elif level > 0:
                block += line

    def parse_f_m_ct(self, blck):

        block_dict = {
            'mae_prop': [],
            'atoms': [],
            'bonds': [],
        }
        current_section = 'mae_prop'
        block = blck.split('\n')
        for line in block:
            if 'm_atom[' in line:
                current_section = 'atoms'
            elif 'm_bond[' in line:
                current_section = 'bonds'
            block_dict[current_section].append(line.strip('\n'))

        return block_dict['atoms'], block_dict['bonds']

    def parse_atoms_block(self, blck):
        res = {}
        separators = []
        cnt = 0
        for line in blck:
            if ':::' in line:
                separators.append(cnt)
            cnt += 1
        columns = blck[2: separators[0]]
        atoms = blck[separators[0]+1:separators[1]]
        x = None
        y = None
        z = None
        Z = None
        for i, c in enumerate(columns):
            c = c.strip(' ')
            if c == 'r_m_x_coord':
                x = i
            elif c == 'r_m_y_coord':
                y = i
            elif c == 'r_m_z_coord':
                z = i
            elif c== 'i_m_atomic_number':
                Z = i
        for line in atoms:
            line = line.replace('"', '').replace('    ', '   ').replace('   ', ' ').replace('  ', ' ').split(' ')[1:]
            atom_Z = int(line[Z+1])
            atom_x = float(line[x+1])
            atom_y = float(line[y+1])
            atom_z = float(line[z+1])
            atom_index = int(line[0])
            res[atom_index] = (atom_Z, atom_x, atom_y, atom_z)

        return res

    def parse_bond_block(self, blck):
        res = {}
        separators = []
        cnt = 0
        for line in blck:
            if ':::' in line:
                separators.append(cnt)
            cnt += 1
        columns = blck[2: separators[0]]
        bonds = blck[separators[0]+1:separators[1]]
        ai = None
        aj = None
        bo = None
        for i, c in enumerate(columns):
            c = c.strip(' ')
            if c == 'i_m_from':
                ai = i
            elif c == 'i_m_to':
                aj = i
            elif c == 'i_m_order':
                bo = i
        for line in bonds:
            line = line.split(' ')[2:]
            atom_i = int(line[ai+1])
            atom_j = int(line[aj+1])
            bond_order = int(line[bo+1])
            res[(atom_i, atom_j)] = bond_order

        return res

    def construct_mol(self, atoms, bonds):
        abi_atoms = {}
        for index, props in atoms.iteritems():
            a = Atom(z=props[0])
            a.x = props[1]
            a.y = props[2]
            a.z = props[3]
            abi_atoms[index] = a
        abi_bonds = []
        for index, bond_order in bonds.iteritems():
            b = Bond(order=bond_order,
                     atom1=abi_atoms[index[0]],
                     atom2=abi_atoms[index[1]]
                     )
            abi_bonds.append(b)
            mol = Molecule()
        mol.add_atoms(abi_atoms.values())
        mol.add_bonds(abi_bonds)
        return mol

    def encode(self, molecules, path):
        output_file = open(path, 'w')
        output_file.write(self.header)

        if isinstance(molecules, Molecule):
            fmct = self.build_fmct(molecules)
            output_file.write(fmct)
        else:
            for molecule in molecules:
                fmct = self.build_fmct(molecule)
                output_file.write(fmct)
        output_file.close()


    def build_fmct(self, molecule):
        fmct_template = open('templates/mae_fmct', 'r').read()
        atom_row_template = '{} 3 {} {} {} 900 2 {} A0A0A0  <>\n'
        bond_row_template = '{} {} {} {}\n'
        natoms = str(len(molecule.atoms))
        nbonds = str(len(molecule.bonds))
        atom_rows_block = []
        cnt = 0
        atom_indexes = {}
        for atom in molecule.atoms:
            cnt += 1
            atom_indexes[atom] = cnt
            s = atom_row_template.format(cnt, atom.x, atom.y, atom.z, atom.Z)
            atom_rows_block.append(s)

        bond_rows_block = []
        cnt = 0
        for bond in molecule.bonds:
            cnt += 1
            atom1, atom2 = (atom for atom in bond)
            index1 = atom_indexes[atom1]
            index2 = atom_indexes[atom2]

            s = bond_row_template.format(cnt, index1, index2, bond.order)
            bond_rows_block.append(s)

        return fmct_template.format(natoms, ''.join(atom_rows_block), nbonds, ''.join(bond_rows_block))

