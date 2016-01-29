# -*- coding:utf-8 -*-
"""

"""
import itertools
import math

from calc.quaternions.quaternion import Quaternion

def translate_mol_to_zero(mol, center):
    x = -1 * center[0]
    y = -1 * center[1]
    z = -1 * center[2]
    for atom in mol.atoms:
        atom.x -= x
        atom.y -= y
        atom.z -= z


def mol_center(mol):
    """
    return center of molecular diameter
    :param mol: molecule
    :return:
    """
    a1, a2, dist = molecule_diameter(mol)
    x = 0.5*(a1.x + a2.x)
    y = 0.5*(a1.y + a2.y)
    z = 0.5*(a1.z + a2.z)
    return x, y, z


def molecule_diameter(mol):
    """
    return longest distance between atoms in molecule (diameter)
    :param mol: molecule
    :return:
    """
    max = -3
    a1max = None
    a2max = None
    for a1, a2 in itertools.combinations(mol.atoms, 2):
        d = distance(a1, a2)
        if d > max:
            max = d
            a1max = a1
            a2max = a2
    return a1max, a2max, max



def distance(a1, a2):
    """
    return distance beetween 2 atoms
    :param a1: atom1
    :param a2: atom2
    :return: float
    """
    x = (a1.x - a2.x)**2
    y = (a1.y - a2.y)**2
    z = (a1.z - a2.z)**2
    return math.sqrt(x+y+z)


def rotate(axis, angle, molecule):
    x, y, z = axis
    qq = Quaternion(0, x, y, z)
    qq.normalize(0.0001)
    theta = angle/2
    w = math.cos(theta)
    x = qq.b * math.sin(theta)
    y = qq.c * math.sin(theta)
    z = qq.d * math.sin(theta)
    q1 = Quaternion(w, x, y, z)
    qc = q1.conjudate()
    n = 0
    for atom in molecule.atoms:
        n += 1
        v = Quaternion(0, atom.x, atom.y, atom.z)
        result_quaternion = q1.mult(v).mult(qc)
        x, y, z = result_quaternion.b, result_quaternion.c, result_quaternion.d
