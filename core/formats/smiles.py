# -*- coding: utf-8 -*-
"""
Specs,Articles :
    http://www.iocd.unam.mx/organica/seminario/weininger88.pdf
    http://www.daylight.com/dayhtml/doc/theory/theory.smiles.html

hydrogen-supressed graph
hydrogen-complete graph
atom : '[' <mass> symbol <chiral> <hcount> <sign<charge>> ']'

"""
import re

from core.exception.exception import BadFormatError
from parser import Parser
from structure.Molecule import Molecule
from structure.Atom import Atom
from structure.Bond import Bond
from utils.periodic_table import periodic_table_by_symbol

SMILES_STRUCTURE = re.compile('^[A-Za-z0-9\[\]\-\=\#\:]*$')
BONDS_CLOSURE = re.compile('^[0-9%]+')
OPEN_BRACKETS = re.compile('^\(+')


class SmilesParser(Parser):
    def __init__(self):
        self.RE_STRUCTURE = re.compile('^(?P<atom>C|\[[a-zA-Z0-9\-+@]*\]|C|N|O|Cl|Br|F|I|S|P|B|\*|n|o|c|s|p)(?P<bonds>[=#%/()0-9\\\\\.]*)$')
        self.SQUARE_BRACKET = re.compile('\[(?P<mass>[0-9]{0,2})(?P<symbol>Th|Rh|[a-gi-zA-GI-Z]{1,2}|H|Hf|Ho|Hg)(?P<chiralsign>@{0,2})(?P<chiralclass>(AL|TB|SP|OH)?[0-9]{0,2})(?P<hpresence>h?|H?)(?P<hcount>[0-9]?)(?P<charge>\+?[0-9]*|-?[0-9]*)\]')
        self.stereo_symbols = {'\\': '+', '/': '-'}

    def is_smiles(self, string):
        """
        Check if it string can be SMILES String or not
        :param string:
        :return:
        """
        if not re.match(SMILES_STRUCTURE, string):
            raise BadFormatError
        return True

    def decode(self, string=''):
        smiles_string = string
        molecule = Molecule()
        previuos_atom = None
        previous_bond = None
        anchor_atom_dct = {}
        atom_stack = []
        after_branch_close = False
        after_branch_close_open = False

        while len(smiles_string) > 0:
            index = 0
            while index < len(smiles_string):
                is_atom = self.RE_STRUCTURE.match(smiles_string[0:index+1])
                if is_atom:
                    atom, smiles_string = smiles_string[0:index+1], smiles_string[index+1:]
                    while smiles_string and self.RE_STRUCTURE.match(atom + smiles_string[0]):
                        atom, smiles_string = atom + smiles_string[0], smiles_string[1:]
                    structure = re.match(self.RE_STRUCTURE, atom).groupdict()
                    atom, bond_expression = structure['atom'], structure['bonds']
                    if atom.startswith('['):
                        # complex atom
                        bracket_structure = self.SQUARE_BRACKET.match(atom).groupdict()
                        atom = bracket_structure['symbol']
                    atom = atom.capitalize()
                    try:
                        atom = Atom(periodic_table_by_symbol[atom]['n'])
                    except Exception as err:
                        pass
                    molecule.atoms.add(atom)
                    if previuos_atom:
                        previous_bond.atoms.add(atom)
                        if after_branch_close:
                            previous_bond.atoms.add(atom_stack.pop())
                            after_branch_close = False
                        if after_branch_close_open:
                            previous_bond.atoms.add(atom_stack[-1])
                            after_branch_close_open = False
                        else:
                            previous_bond.atoms.add(previuos_atom)
                        molecule.bonds.add(previous_bond)
                    previuos_atom = atom
                    if not bond_expression:
                        bond = Bond(order=1)
                    else:
                        # we have something in bond
                        # we want splt it on 2 part first part bond with number, other all except it
                        S = re.compile('^(?P<numbering>.*\d+)(?P<tail>.*)$')
                        s = S.match(bond_expression)
                        if s:
                            s = s.groupdict()
                            numbering, tail = s['numbering'], s['tail']
                            if not numbering[0].isdigit():
                                cis_trans_sign, numbers = numbering[0], numbering[1:]
                            else:
                                cis_trans_sign, numbers = None, numbering
                            numbers_list = []
                            while numbers:
                                if not numbers[0] == '%':
                                    number, numbers = numbers[0], numbers[1:]
                                else:
                                    number, numbers = numbers[1:3], numbers[3:]
                                numbers_list.append(number)
                            # DO SOMETHING WITH THAT BOND
                        else:
                            tail = bond_expression

                        if not tail:
                            # tail = brancing + bond
                            bond = Bond(order=1)
                        else:
                            T = re.compile('^(?P<branch>[()]{0,2})(?P<cis_trans_sign>[/\\\\]?)(?P<bond_type>.?)$')
                            t = T.match(tail).groupdict()
                            branch, cis_trans_sign, bond_type = t['branch'], t['cis_trans_sign'], t['bond_type']
                            if branch:
                                if branch == '(':
                                    atom_stack.append(atom)
                                elif branch == ')':
                                    after_branch_close = True
                                elif branch == ')(':
                                    after_branch_close_open = True
                                else:
                                    print branch
                            cis_trans_sign = cis_trans_sign or None
                            if not bond_type:
                                bond = Bond(order=1)
                            elif bond_type == '=':
                                bond = Bond(order=2)
                            elif bond_type == '#':
                                bond = Bond(order=3)
                            elif bond_type == '.':
                                # ionic character of bond
                                bond = Bond(order=1)
                            else:
                                print bond_type

                                bond = Bond(order=1)
                    previous_bond = bond

                else:
                    index += 1
        return molecule


    def encode(self, molecule):
        return None


    def get_canonical(self):
        return self.encode(self.decode(''))


if __name__ == '__main__':
    p = SmilesParser()
    import os
    import time
    n = 0
    s = time.time()
    smiles = open(os.getcwd() + '/smiles.txt', 'r')
    for smile in smiles:
        mol = p.decode(smile.split(' ')[0])
        n += 1
        if n % 10000 == 0:
            print n, (time.time() - s)/n

