# -*- coding: utf-8 -*-
"""
Specs,Articles :
    http://www.iocd.unam.mx/organica/seminario/weininger88.pdf
    http://www.daylight.com/dayhtml/doc/theory/theory.smiles.html

hydrogen-supressed graph
hydrogen-complete graph
atom : '[' <mass> symbol <chiral> <hcount> <sign<charge>> ']'

parser.decode(SMILES_STRING) --> molecule internal representation
# TODO! aromacity order of bond
# add more simplicity

"""
import re

from core.exception.exception import BadFormatError
from parser import Parser
from structure.Molecule import Molecule
from structure.Atom import Atom
from structure.Bond import Bond
from utils.periodic_table import periodic_table_by_symbol, periodic_table


class SmilesParser(Parser):
    def __init__(self):
        self.SMILES_STING = re.compile('^[A-Za-z0-9\[\]\-=#:\\\\/().\+@%\*]*$')
        self.RE_STRUCTURE = re.compile('(?P<atom>C|\[[a-zA-Z0-9\-+@]*\]'
                                       '|C|N|O|Cl|Br|F|I|S|P|B|\*|n|o|c|s|p)'
                                       '(?P<bonds>[=#%/()0-9\\\\\.]*)')
        self.SQUARE_BRACKET = re.compile('\[(?P<mass>[0-9]{0,3})'
                                         '(?P<symbol>Th|Rh|[a-gi-zA-GI-Z]{1,2}|H|Hf|Ho|Hg)'
                                         '(?P<chiralsign>@{0,2})'
                                         '(?P<chiralclass>(AL|TB|SP|OH)?[0-9]{0,2})'
                                         '(?P<hpresence>h?|H?)(?P<hcount>[0-9]?)'
                                         '(?P<charge>\+?[0-9]*|-?[0-9]*)\]')
        self.BREAKING_BOND = re.compile('^(?P<breaking_bond>.*\d+)(?P<tail>.*)$')
        self.TAIL = re.compile('^(?P<branch>[()]{0,2})(?P<cis_trans_sign>[/\\\\]?)(?P<bond_type>.?)$')

    def is_smiles(self, string):
        """
        Check if it string can be SMILES String or not
        :param string:
        :return: bool
        """
        if not self.SMILES_STING.match(string):
            raise BadFormatError
        return True

    def decode(self, string):
        self.is_smiles(string)
        molecule = Molecule()
        previous_atom = None
        previous_bond = None
        atom_stack = []
        numbering_stack = {}
        after_branch_close = False
        after_branch_close_open = False

        parsing_string = re.findall(self.RE_STRUCTURE, string)
        for atom_exp, bond_exp in parsing_string:
            #PARSING ATOM EXPRESSION
            if atom_exp[0] != '[':
                atom_symbol = atom_exp
                neuthron_number = 0
                charge = 0
                h_count = 0
                chiral_sign = ''
                chiral_class = ''
                bracket = False
            else:
                complex_atom = self.SQUARE_BRACKET.match(atom_exp).groupdict()
                atom_symbol = complex_atom['symbol']
                neuthron_number = int(complex_atom['mass'] or 0)
                h_count = 0
                if complex_atom['hpresence']:
                    h_count = int(complex_atom['hcount'] or 1)

                charge_string = complex_atom['charge']
                charge = 0
                if charge_string:
                    if len(charge_string) == 1:
                        if charge_string == '-':
                            charge = -1
                        elif charge_string == '+':
                            charge = 1
                    else:
                        if charge_string[0] == '-':
                            charge = 0 - int(charge_string[1:])
                        else:
                            charge = int(charge_string[1:])
                chiral_sign = complex_atom['chiralsign']
                chiral_class = complex_atom['chiralclass']
                bracket = True

            if atom_symbol.islower():
                aromatic = True
            else:
                aromatic = False
            atom_symbol = atom_symbol.capitalize()
            z = periodic_table_by_symbol[atom_symbol]['z']

            # CREATING ATOM
            atom = Atom(z=z, n=neuthron_number or z, charge=charge, chiral=chiral_class + chiral_sign, aromatic=aromatic, bracket=bracket)
            molecule.atoms.append(atom)
            for x in range(h_count):
                h = Atom(1)
                molecule.atoms.append(h)
                bond_h = Bond(order=1, atom1=atom, atom2=h)
                molecule.bonds.append(bond_h)

            if previous_atom:
                previous_bond.add_atom(atom)
                if after_branch_close:
                    atom_pop = atom_stack.pop()
                    previous_bond.add_atom(atom_pop)
                    after_branch_close = False
                if after_branch_close_open:
                    previous_bond.add_atom(atom_stack[-1])
                    after_branch_close_open = False
                else:
                    previous_bond.add_atom(previous_atom)
                molecule.bonds.append(previous_bond)
            previous_atom = atom
            if not bond_exp:
                bond = Bond(order=1)
            else:
                split_bond_expression = self.BREAKING_BOND.match(bond_exp)
                if split_bond_expression:
                    split_bond_expression = split_bond_expression.groupdict()
                    breaking_bond = split_bond_expression['breaking_bond']
                    tail = split_bond_expression['tail']
                    if not breaking_bond[0].isdigit() and breaking_bond[0] != '%':
                        cis_trans_sign, numbers = breaking_bond[0], breaking_bond[1:]
                    else:
                        cis_trans_sign, numbers = None, breaking_bond
                    numbers_list = []
                    while numbers:
                        if not numbers[0] == '%':
                            number, numbers = numbers[0], numbers[1:]
                        else:
                            number, numbers = numbers[1:3], numbers[3:]
                        numbers_list.append(number)

                    for number in numbers_list:
                        if not number in numbering_stack:
                            numbering_stack[number] = atom
                        else:
                            numbering_bond = Bond(order=1, cis_trans=cis_trans_sign, atom1=numbering_stack[number], atom2=atom)
                            molecule.bonds.append(numbering_bond)
                            del(numbering_stack[number])
                else:
                    tail = bond_exp
                if not tail:
                    bond = Bond(order=1)
                else:
                    t = self.TAIL.match(tail).groupdict()
                    branch, cis_trans_sign, bond_type = t['branch'], t['cis_trans_sign'], t['bond_type']
                    if branch:
                        if branch == '(':
                            atom_stack.append(atom)
                        elif branch == ')':
                            after_branch_close = True
                        elif branch == ')(':
                            after_branch_close_open = True
                        cis_trans_sign = cis_trans_sign or None
                    if not bond_type:
                        bond = Bond(order=1, cis_trans=cis_trans_sign)
                    elif bond_type == '=':
                        bond = Bond(order=2, cis_trans=cis_trans_sign)
                    elif bond_type == '#':
                        bond = Bond(order=3, cis_trans=cis_trans_sign)
                    elif bond_type == '.':
                        # ionic character of bond
                        bond = Bond(order=1, cis_trans=cis_trans_sign)

            previous_bond = bond


        # add atoms h on unocuppy atoms
        free_atoms = {}
        for index, atom in enumerate(molecule.atoms):
            #print atom
            if atom.bracket:
                continue
            aromatic = 1 if atom.aromatic else 0
            if periodic_table[atom.Z]['symbol'] == 'C':
                free_bond = 4 - sum([bond.order for bond in atom.bonds])
                if free_bond > 0:
                    free_bond -= aromatic
                    free_atoms[atom] = free_bond
            elif periodic_table[atom.Z]['symbol'] == 'O':
                free_bond = 2 - sum([bond.order for bond in atom.bonds])
                if free_bond > 0:
                    free_atoms[atom] = free_bond
            elif periodic_table[atom.Z]['symbol'] == 'N':
                if sum([bond.order for bond in atom.bonds]) <= 3:
                    free_bond = 3 - sum([bond.order for bond in atom.bonds]) - atom.charge - aromatic
                else:
                    free_bond = 5 - sum([bond.order for bond in atom.bonds]) - atom.charge
                if free_bond:
                    free_atoms[atom] = free_bond
            elif periodic_table[atom.Z]['symbol'] == 'P':
                if sum([bond.order for bond in atom.bonds]) <= 3:
                    free_bond = 3 - sum([bond.order for bond in atom.bonds]) - atom.charge
                else:
                    free_bond = 5 - sum([bond.order for bond in atom.bonds]) - atom.charge
                if free_bond:
                    free_atoms[atom] = free_bond
            elif periodic_table[atom.Z]['symbol'] == 'S':
                if sum([bond.order for bond in atom.bonds]) <= 2:
                    free_bond = 2 - sum([bond.order for bond in atom.bonds]) - atom.charge
                elif sum([bond.order for bond in atom.bonds]) <= 4:
                    free_bond = 4 - sum([bond.order for bond in atom.bonds]) - atom.charge
                else:
                    free_bond = 6 - sum([bond.order for bond in atom.bonds]) - atom.charge
                if free_bond:
                    free_atoms[atom] = free_bond

        for atom, h_number in free_atoms.iteritems():
            for x in range(h_number):
                h = Atom(1)
                molecule.atoms.append(h)
                b = Bond(order=1, atom1=atom, atom2=h)
                molecule.bonds.append(b)
                b.add_atom(h)
                b.add_atom(atom)






        return molecule


    def encode(self, molecule):
        return None


    def get_canonical(self):
        return self.encode(self.decode(''))


if __name__ == '__main__':
    import os
    p = SmilesParser()
    smiles = open(os.getcwd() + '/smiles.txt', 'r')
    for smile in smiles:
        mol = p.decode(smile.split(' ')[0])




















