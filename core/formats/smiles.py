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

from core.exception.exception import BadFormatError, ParseError
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
        self.RE_STRUCTURE = '^(?P<atom>\[[a-zA-Z0-9\-+@]*\]|C|N|O|Cl|Br|F|I|S|P|B|\*|n|o|c|s|p)(?P<stereo1>/?|\\\\?)(?P<bond1>|=?|\.|#?|/[0-9]*)(?P<number1>[0-9%]*)(?P<branching>[()]*)(?P<stereo2>/?|\\\\?)(?P<number2>[0-9%]*)(?P<bond2>/?|=?|\\\\?|\.|#?|/[0-9]*)$'
        self.SQUARE_BRACKET = '\[(?P<mass>[0-9]{0,2})(?P<symbol>Th|Rh|[a-gi-zA-GI-Z]{1,2}|H|Hf|Ho|Hg)(?P<chiralsign>@{0,2})(?P<chiralclass>(AL|TB|SP|OH)?[0-9]{0,2})(?P<hpresence>h?|H?)(?P<hcount>[0-9]?)(?P<charge>\+?[0-9]*|-?[0-9]*)\]'

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
        #self.is_smiles(string)
        def test_pattern(pattern):
            # Structure (atom,stereo1,bond1,number1,branching,stereo2,bond2,number2)
            if re.match(self.RE_STRUCTURE, pattern):
                return True
            return False


        smiles_string = string
        atom_stack = []
        molecule = Molecule()

        while len(smiles_string) > 0:
            index = 0
            while index < len(smiles_string):
                is_atom = test_pattern(smiles_string[0:index+1])
                if is_atom:
                    atom, smiles_string = smiles_string[0:index+1], smiles_string[index+1:]
                    while smiles_string and test_pattern(atom + smiles_string[0]):
                        atom, smiles_string = atom + smiles_string[0], smiles_string[1:]
                    structure = re.match(self.RE_STRUCTURE, atom)
                    structure = structure.groupdict()
                    atom = structure['atom']
                    if atom.startswith('['):
                        atom_in_bracket = re.match(self.SQUARE_BRACKET, atom).groupdict()
                        if not atom_in_bracket:
                            print atom
                        atom = atom_in_bracket['symbol']
                    try:
                        atom = atom.capitalize()
                        atom = Atom(periodic_table_by_symbol[atom]['n'])
                        molecule.atoms.add(atom)
                    except Exception as err:
                        pass
                    atom_stack.append(atom)
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
        l = p.decode(smile.split(' ')[0])
























