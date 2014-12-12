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
from utils.periodic_table import periodic_table
periodic_dct = {}

for k, v in periodic_table.items():
    periodic_dct[v['symbol']] = k

SMILES_STRUCTURE = re.compile('^[A-Za-z0-9\[\]\-\=\#\:]*$')
SQUARE_BRACKET = re.compile('\[(?P<mass>[0-9]{0,2})(?P<symbol>[a-zA-GI-Z]{1,2}|H|Hf|Ho|Hg)(?P<chiralsign>@{0,2})(?P<chiralclass>(AL|TB|SP|OH)?[0-9]{0,2})(?P<hcount>H?[0-9]?)(?P<charge>\+?[0-9]*|-?[0-9]*)\]')
BONDS_CLOSURE = re.compile('^[0-9%]+')
OPEN_BRACKETS = re.compile('^\(+')


class SmilesParser(Parser):
    def __init__(self):
        pass

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
            RE_STRUCTURE = '^(?P<atom>\[[a-zA-Z0-9\-+@]*\]|O|C|N|Cl|Br|F|I|S|P|B|\*|n|o|c|s|p)(?P<stereo1>/?|\\\\?)(?P<bond1>|=?|\.|#?|/[0-9]*)(?P<number1>[0-9%]*)(?P<branching>[()]*)(?P<stereo2>/?|\\\\?)(?P<number2>[0-9%]*)(?P<bond2>/?|=?|\\\\?|\.|#?|/[0-9]*)$'
            # Structure (atom,stereo1,bond1,number1,branching,stereo2,bond2,number2)
            if re.match(RE_STRUCTURE, pattern):
                return True
            return False


        smiles_string = string
        structure_stack = []

        while len(smiles_string) > 0:
            index = 0
            atom = ''
            while index < len(smiles_string):
                is_atom = test_pattern(smiles_string[0:index+1])
                if is_atom:
                    atom, smiles_string = smiles_string[0:index+1], smiles_string[index+1:]
                    while smiles_string and test_pattern(atom + smiles_string[0]):
                        atom, smiles_string = atom + smiles_string[0], smiles_string[1:]
                    structure_stack.append(atom)
                else:
                    index += 1
        #print structure_stack
        return None


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


























