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


SMILES_STRUCTURE = re.compile('^[A-Za-z0-9\[\]\-\=\#\:]*$')


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
        print(string)
        return None


    def encode(self, molecule):
        return None


    def get_canonical(self):
        return self.encode(self.decode(''))


if __name__ == '__main__':
    p = SmilesParser()
