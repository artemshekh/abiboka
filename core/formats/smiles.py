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
                        previous_bond.atoms.add(previuos_atom)
                        molecule.bonds.add(previous_bond)
                    previuos_atom = atom
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
    from collections import defaultdict
    d = defaultdict(lambda: 0)
    smiles = open(os.getcwd() + '/smiles.txt', 'r')
    for smile in smiles:
        mol = p.decode(smile.split(' ')[0])
