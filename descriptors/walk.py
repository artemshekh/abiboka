# -*- coding: utf-8 -*-
"""
Walk and path counts
"""

from descriptor_utils import path_sequence_matrix

def mpc(molecule, order=1):
    """
    molecular path order
    :param molecule:
    :param order:
    :return:
    """
    m = path_sequence_matrix(molecule, order)
    return sum([row[-1] for row in m.matrix])/2.0

def tpc(molecule):
    m = path_sequence_matrix(molecule)
    return len(molecule.atoms) + sum([sum(row) for row in m.matrix])/2.0
