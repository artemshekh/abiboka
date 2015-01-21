# -*- coding: utf-8 -*-
"""
InChi parser
IUPAC International Chemical Identifier


"""
from parser import Parser
import re

from structure.Molecule import Molecule
from structure.Atom import Atom
from utils.periodic_table import periodic_table_by_symbol

RE_HILL = re.compile(('(?P<atom>[A-Z][a-z]?)(?P<number>[0-9]*)'))

class InChIParser(Parser):

    def decode(self, string=''):
        print string
        molecule = Molecule()
        sections = string.split('/')
        hill_formula = sections[1]
        #print hill_formula
        atoms = re.findall(RE_HILL, hill_formula)

        for atom in atoms:
            c = int(atom[1]) if atom[1] else 1
            for count in range(c):
                a = Atom(z=periodic_table_by_symbol[atom[0]]['z'])
                molecule.atoms.append(a)

        #print molecule.atoms, molecule.hydrogen_suppressed.atoms


        return None