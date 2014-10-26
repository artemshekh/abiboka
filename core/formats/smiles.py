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
        self.is_smiles(string)
        #new_molecule = Molecule()
        #Find root atom remove it from string with all parameters
        #find bond find next atom
        #

        print(string)
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
    n = 0

    def parse_atom(smile):
        atomstring = ''
        if smile[0:2] in ['Br', 'Cl']:
            atomstring = atomstring + smile[0:2]
            smile = smile[2:]
            #print(atom, smile)
            #atom = Atom(z=periodic_dct[atom])
        elif smile[0] in ['C', 'N', 'O', 'P', 'S', 'F', 'I', 'B', 'c', 'n', 's', 'o', '*', 'p']:
            atomstring = atomstring + smile[0]
            symbol = smile[0]
            smile = smile[1:]
            #print(atom, smile)
            #atom = Atom(z=periodic_dct[symbol.capitalize()])
        elif smile[0] == '[':
            atom = re.search(SQUARE_BRACKET, smile).groupdict()
            atomstring = atomstring + smile[:re.search(SQUARE_BRACKET, smile).end()]
            smile = smile[re.search(SQUARE_BRACKET, smile).end():]
            #atom = Atom(z=periodic_dct[atom['symbol'].capitalize()])
            #print(atom, smile)

        else:
            raise ParseError
        return atomstring, smile

    def parse_bond(string):
        bondstring = ''
        if string[0] == '/':
            bondstring += string[0]
            string = string[1:]

        if string[0] == '\\':
            bondstring += string[0]
            string = string[1:]
            if string[0] == '\\':
                bondstring += string[0]
                string = string[1:]
        if re.search(BONDS_CLOSURE, string):
            bondstring += string[:re.search(BONDS_CLOSURE, string).end()]
            string = string[re.search(BONDS_CLOSURE, string).end():]
        if not string:
            return bondstring, string
        if string[0] == ')':
                bondstring += string[0]
                branch_stack.pop()
                string = string[1:]
        if string[0] == '(':
                print()
                bondstring += string[0]
                branch_stack.append(atom)
                string = string[1:]

        if string[0] == '/':
            bondstring += string[0]
            string = string[1:]
        if string[0] == '\\':
            bondstring += string[0]
            string = string[1:]
            if string[0] == '\\':
                bondstring += string[0]
                string = string[1:]
        if string[0] == '=':
            bondstring += string[0]
            bond = Bond(order=2)
            string = string[1:]
        elif string[0] == '#':
            bondstring += string[0]
            bond = Bond(order=3)
            string = string[1:]
        else:
            bond = Bond(order=1)
        if string[0] == '.':
            bondstring += string[0]
            string = string[1:]
        if re.search(BONDS_CLOSURE, string):
            bondstring += string[:re.search(BONDS_CLOSURE, string).end()]
            string = string[re.search(BONDS_CLOSURE, string).end():]
        if string[0] == ')':
                bondstring += string[0]
                branch_stack.pop()
                string = string[1:]
        if string[0] == '(':
                print()
                bondstring += string[0]
                branch_stack.append(atom)
                string = string[1:]

        if len(bondstring) >1:
            print(bondstring)






        return bondstring, string
    n = 0
    for line in smiles:
        n +=1
        if n <6454701:
            continue
        atoms = []
        branch_stack = []
        smile = line.split(' ')[0]
        while smile:
            print(smile)
            atom, smile = parse_atom(smile)
            atoms.append(atom)
            #print(atom, smile, line)
            #print(type(atom))
            if smile:
                if smile[0] == ')':
                    branch_stack.pop()
                    smile = smile[1:]
                if smile[0] == '(':
                    branch_stack.append(atom)
                    smile = smile[1:]

                bond, smile = parse_bond(smile)
            import string as s_
            #if smile and smile[0] not in s_.ascii_letters + '[*':
            #    d[smile[0:20]] += 1
            #print atom
            print(atom, bond)
        print(n)





    d = d.items()
    d.sort(key=lambda x: -x[1])
    for k,v in d:
        print(k,v)






















