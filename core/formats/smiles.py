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

SMILES_STRUCTURE = re.compile('^[A-Za-z0-9\[\]\-\=\#\:]*$')


def is_smiles(string):
    if not re.match(SMILES_STRUCTURE, string):
        raise BadFormatError
    return True


def decode(string=''):
    is_smiles(string)
    print(string)
    return None


def encode(molecule):
    return None


def get_canonical():
    return encode(decode(''))


if __name__ == '__main__':
    decode('[c]=c')
