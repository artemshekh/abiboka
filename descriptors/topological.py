# -*- coding: utf-8 -*-
"""
Topological descriptors
"""
import operator
from collections import Counter
from utils.periodic_table import periodic_table
from calc.matrixes.matrix import AdjacencyMatrix, Matrix
from descriptors.descriptor_utils import path_sequence_matrix, walk_vector
from descriptors.walk import mpc
from descriptors.ring_descriptor import cyclomatic_number
import math


def zm1(molecule):
    return sum([len(atom.bonds)**2 for atom in molecule.hydrogen_suppressed.atoms])

def zm1_H(molecule):
    return sum([len(atom.bonds)**2 for atom in moleculeatoms])

close_shell = [1, 3, 11, 19, 37, 55, 87]

def valence_electrones(atom):
    valence_e = 4
    for period,n in enumerate(close_shell):
        if atom.Z < n:
            core_e, valence_e = close_shell[period-1]-1, atom.Z - close_shell[period-1] + 1
            break
    if valence_e:
       return valence_e
    else:
        return atom.Z - 87


def valence_degree(atom):
    # all valence electrons of the ith atom
    vd = valence_electrones(atom)
    for _ in atom.connected_with():
        if atom.Z == 1:
            vd -= 1
    return vd

def valence_degree_(atom):
    # the electronic identity of the atom in terms of both valence electron and core electron counts
    vd = valence_electrones(atom)
    for _ in atom.connected_with():
        if atom.Z == 1:
            vd -= 1
    vd = float(vd)/(atom.Z - valence_e - 1)
    return vd


def zm1v(molecule):
    return sum([valence_degree(atom)**2 for atom in molecule.hydrogen_suppressed.atoms])

def zm1v_(molecule):
    # depends from core electrons
    return sum([valence_degree(atom)**2 for atom in molecule.hydrogen_suppressed.atoms])

def kupchik_degree(atom):
    return  (float(periodic_table[6]['covalent_radius'])/periodic_table[atom.Z]["covalent_radius"])* valence_degree(atom)

def zm1kup(molecule):
    return sum([kupchik_degree(atom)**2 for atom in molecule.atoms])

def madan_degree(atom):
    return sum([periodic_table[atom1.Z]["relative_atomic_weight"]/periodic_table[6]["relative_atomic_weight"] for atom1 in atom.connected_with()])

def zm1mad(molecule):
    return sum([madan_degree(atom) ** 2 for atom in molecule.atoms])

def permutation_additive(atom, perm_coefficient):
    return valence_degree(atom) + sum([perm_coefficient * valence_degree(atom) for _ in atom.connected_with()])

def permutation_multiplicative(atom, perm_coefficient):
    return valence_degree(atom) + reduce(operator.mul, [perm_coefficient * valence_degree(atom) for _ in atom.connected_with()])

def zm1per(molecule, permutation_coefficient):
    return sum([permutation_additive(atom, permutation_coefficient) ** 2 for atom in molecule.atoms])

def zm1mulper(molecule, permutation_coefficient):
    return sum([permutation_multiplicative(atom, permutation_coefficient) ** 2 for atom in molecule.atoms])

def zm2(molecule):
    _ = []
    for bond in molecule.hydrogen_suppressed.bonds:
        atoms = [atom for atom in bond]
        _.append(len(atoms[0].bonds)*len(atoms[1].bonds))
    return sum(_)

def zm2v(molecule):
    _ = []
    for bond in molecule.bonds:
        atoms = [atom for atom in bond]
        _.append(valence_degree(atoms[0])*valence_degree(atoms[1]))
    return sum(_)

def zm2kup(molecule):
    _ = []
    for bond in molecule.bonds:
        atoms = [atom for atom in bond]
        _.append(kupchik_degree(atoms[0])*kupchik_degree(atoms[1]))
    return sum(_)


def zm2mad(molecule):
    _ = []
    for bond in molecule.bonds:
        atoms = [atom for atom in bond]
        _.append(madan_degree(atoms[0])*madan_degree(atoms[1]))
    return sum(_)

def zm2per(molecule, permutative_coefficent):
    _ = []
    for bond in molecule.bonds:
        atoms = [atom for atom in bond]
        _.append(permutation_additive(atoms[0], permutative_coefficent)*permutation_additive(atoms[1], permutative_coefficent))
    return sum(_)

def zm2mulper(molecule, permutative_coefficent):
    _ = []
    for bond in molecule.bonds:
        atoms = [atom for atom in bond]
        _.append(permutation_multiplicative(atoms[0], permutative_coefficent)*permutation_multiplicative(atoms[1], permutative_coefficent))
    return sum(_)

def on0(molecule):
    return sum([1.0/len(atom.bonds) for atom in molecule.hydrogen_suppressed().atoms])

def on0v(molecule):
    _ = []
    for atom in molecule.atoms:
        v = valence_degree(atom)
        if v:
            _.append(1.0/v)
    return sum(_)

def on1(molecule):
    _ = []
    for bond in molecule.hydrogen_suppressed.bonds:
        atoms = [atom for atom in bond]
        _.append(1.0/len(atoms[0].bonds)*1.0/len(atoms[1].bonds))
    return sum(_)

def on1v(molecule):
    _ = []
    for bond in molecule.hydrogen_suppressed.bonds:
        atoms = [atom for atom in bond]
        _.append(1.0/valence_degree(atoms[0])*1.0/valence_degree(atoms[1]))
    return sum(_)

def qindex(molecule):
    molecule = molecule.hydrogen_suppressed
    return 3 -  2* len(molecule.atoms) + zm1(molecule)/2.0

def bbi(molecule):
    return sum([len(atom.bonds)*(len(atom.bonds) - 1) for atom in molecule.hydrogen_suppressed.atoms])/2.0

def snar(molecule):
    return reduce(operator.mul, [len(atom.bonds) for atom in molecule.hydrogen_suppressed.atoms])

def hnar(molecule):
    molecule = molecule.hydrogen_suppressed
    return len(molecule.atoms)/sum([1.0/len(atom.bonds) for atom in molecule.hydrogen_suppressed.atoms])

def gnar(molecule):
    molecule = molecule.hydrogen_suppressed
    return math.pow(snar(molecule), 1.0/len(molecule.atoms))

def xt(molecule):
    molecule = molecule.hydrogen_suppressed
    return 1.0/math.sqrt(snar(molecule))

def dz(molecule):
    _ = []
    p = 2
    for atom in molecule.atoms:
        for period,n in enumerate(close_shell):
            if atom.Z < n:
                p = period
                break
        _.append(valence_degree(atom)/float(p))
    return sum(_)

def ram(molecule):
    s = 0
    for atom in molecule.hydrogen_suppressed.atoms:
        if len(atom.bonds) > 2:
            s += (len(atom.bonds) -2)
    return s

def bli(molecule):
    _ = []
    molecule = molecule.hydrogen_suppressed
    for bond in molecule.bonds:
        atoms = [atom for atom in bond]
        _.append(1/math.sqrt(valence_degree(atoms[0])*valence_degree(atoms[1])))
    return sum(_)/6

def pol(molecule):
    # Wiener polarity number - polarity index
    molecule = molecule.hydrogen_suppressed
    matrix = AdjacencyMatrix.from_molecule(molecule)
    m = matrix.matrix
    for i, row in enumerate(m):
        for j, value in enumerate(row):
            if i!=j and value == 0:
                m[i][j] = 1000
    for k in range(matrix.rows()):
        for i in range(matrix.rows()):
            for j in range(matrix.rows()):
                m[i][j] = min(m[i][j], m[i][k] + m[k][j])
    p = 0
    for row in m:
        for v in row:
            if v == 3:
                p+=1
    return p/2

def prs(molecule):
    return reduce(operator.mul, [len(atom.bonds) for atom in molecule.atoms])

def lprs(molecule):
    return math.log(prs(molecule))

def msd(molecule):
    molecule = molecule.hydrogen_suppressed
    matrix = AdjacencyMatrix.from_molecule(molecule)
    m = matrix.matrix
    for i, row in enumerate(m):
        for j, value in enumerate(row):
            if i!=j and value == 0:
                m[i][j] = 1000
    for k in range(matrix.rows()):
        for i in range(matrix.rows()):
            for j in range(matrix.rows()):
                m[i][j] = min(m[i][j], m[i][k] + m[k][j])
    _ = 0
    for row in m:
        for v in row:
            _ += v*v
    return math.sqrt(float(_)/(len(molecule.atoms)*(len(molecule.atoms) - 1)))

def spi(molecule):
    molecule = molecule.hydrogen_suppressed
    pendant_vertexes = []
    for index, atom in enumerate(molecule.atoms):
        if len(atom.bonds) == 1:
            pendant_vertexes.append(index)
    if not pendant_vertexes:
        return 0
    matrix = AdjacencyMatrix.from_molecule(molecule)
    m = matrix.matrix
    for i, row in enumerate(m):
        for j, value in enumerate(row):
            if i!=j and value == 0:
                m[i][j] = 1000
    for k in range(matrix.rows()):
        for i in range(matrix.rows()):
            for j in range(matrix.rows()):
                m[i][j] = min(m[i][j], m[i][k] + m[k][j])
    descriptor = 0
    for row in m:
        descriptor += reduce(operator.mul, [row[x] for x in pendant_vertexes])
    descriptor = math.sqrt(descriptor)
    return descriptor

def distance_matrix(molecule):
    matrix = AdjacencyMatrix.from_molecule(molecule)
    m = matrix.matrix
    for i, row in enumerate(m):
        for j, value in enumerate(row):
            if i!=j and value == 0:
                m[i][j] = 1000
    for k in range(matrix.rows()):
        for i in range(matrix.rows()):
            for j in range(matrix.rows()):
                m[i][j] = min(m[i][j], m[i][k] + m[k][j])
    return m


def pji2(molecule):
    molecule = molecule.hydrogen_suppressed
    m = distance_matrix(molecule)
    diametr, radius = 0, 100000
    for row in m:
        _ = max(row)
        if _ > diametr:
            diametr = _
        if _ < radius:
            radius = _
    return (diametr - radius)/float(diametr)

def ecc(molecule):
    molecule = molecule.hydrogen_suppressed
    m = distance_matrix(molecule)
    return sum([max(row) for row in m])

def aecc(molecule):
    molecule = molecule.hydrogen_suppressed
    m = distance_matrix(molecule)
    return sum([max(row) for row in m])/float(len(molecule.atoms))

def decc(molecule):
    molecule = molecule.hydrogen_suppressed
    m = distance_matrix(molecule)
    aecc_ = sum([max(row) for row in m])/float(len(molecule.atoms))
    return sum([abs(max(row) - aecc_) for row in m])/float(len(molecule.atoms))

def agdd(molecule):
    molecule = molecule.hydrogen_suppressed
    m = distance_matrix(molecule)
    return sum([sum(row) for row in m])/float(len(molecule.atoms))

def mddd(molecule):
    molecule = molecule.hydrogen_suppressed
    m = distance_matrix(molecule)
    agdd_ = sum([sum(row) for row in m])/float(len(molecule.atoms))
    return sum([abs(sum(row) - agdd_) for row in m])/float(len(molecule.atoms))

def unip(molecule):
    molecule = molecule.hydrogen_suppressed
    m = distance_matrix(molecule)
    return min([sum(row) for row in m])

def cent(molecule):
    molecule = molecule.hydrogen_suppressed
    m = distance_matrix(molecule)
    unip_ = min([sum(row) for row in m])
    return sum([sum(row) for row in m]) - (len(molecule.atoms) * unip_)

def var(molecule):
    molecule = molecule.hydrogen_suppressed
    m = distance_matrix(molecule)
    unip_ = min([sum(row) for row in m])
    return max([sum(row) - unip_ for row in m])

def icr(molecule):
    molecule = molecule.hydrogen_suppressed
    m = distance_matrix(molecule)
    c = Counter()
    for row in m:
        c[max(row)] += 1
    descriptor = 0
    for k, v in c.iteritems():
        ng = float(v)/len(molecule.atoms)
        descriptor += ng * (math.log(ng, 2))
    return -1 * descriptor

def smti(molecule):
    molecule = molecule.hydrogen_suppressed
    ad = AdjacencyMatrix.from_molecule(molecule) + Matrix(distance_matrix(molecule))
    vertex_degree = Matrix([[x] for x in [len(atom.bonds) for atom in molecule.atoms]])
    descriptor = 0
    _ = ad * vertex_degree
    for row in _.matrix:
        descriptor += sum(row)
    return descriptor

def smtiv(molecule):
    molecule = molecule.hydrogen_suppressed
    ad = AdjacencyMatrix.from_molecule(molecule) + Matrix(distance_matrix(molecule))
    vertex_degree = Matrix([[x] for x in [valence_degree(atom) for atom in molecule.atoms]])
    descriptor = 0
    _ = ad * vertex_degree
    for row in _.matrix:
        descriptor += sum(row)
    return descriptor

def gmti(molecule):
    molecule = molecule.hydrogen_suppressed
    m = distance_matrix(molecule)
    vertex_degree = [len(atom.bonds) for atom in molecule.atoms]
    descriptor = 0
    for i, row in enumerate(m):
        for j, value in enumerate(row):
            descriptor += value * vertex_degree[i]*vertex_degree[j]
    return descriptor

def gmtiv(molecule):
    molecule = molecule.hydrogen_suppressed
    m = distance_matrix(molecule)
    vertex_degree = [valence_degree(atom) for atom in molecule.atoms]
    descriptor = 0
    for i, row in enumerate(m):
        for j, value in enumerate(row):
            descriptor += value * vertex_degree[i]*vertex_degree[j]
    return descriptor

def xu(molecule):
    molecule = molecule.hydrogen_suppressed
    a = math.sqrt(len(molecule.atoms))
    vertex_degree = [len(atom.bonds) for atom in molecule.atoms]
    m = distance_matrix(molecule)
    distance_degree = [sum(row) for row in m]
    numerator, denumerator = 0, 0
    for i, v in enumerate(vertex_degree):
        numerator += v * distance_degree[i]**2
        denumerator += v * distance_degree[i]
    return a *  math.log(float(numerator)/denumerator)

def csi(molecule):
    molecule = molecule.hydrogen_suppressed
    vertex_degree = [len(atom.bonds) for atom in molecule.atoms]
    m = distance_matrix(molecule)
    max_distance = [max(row) for row in m]
    return sum([x[0]*x[1] for x in zip(vertex_degree, max_distance)])

def wap(molecule):
    molecule = molecule.hydrogen_suppressed
    descriptor = 0
    for row in m.matrix:
        for i, value in enumerate(row):
            descriptor += value * (i+1)
    return descriptor/2 + len(molecule.atoms)

def s1k(molecule):
    molecule = molecule.hydrogen_suppressed
    alpha = sum([(float(periodic_table[atom.Z]['covalent_radius'])/periodic_table[6]['covalent_radius']) -1 for atom in molecule.atoms])
    a = len(molecule.atoms)
    p = len(molecule.bonds)/2
    return float(((a + alpha) * ((a + alpha -1)**2)))/((p + alpha)**2)

def s2k(molecule):
    molecule = molecule.hydrogen_suppressed
    alpha = sum([(float(periodic_table[atom.Z]['covalent_radius'])/periodic_table[6]['covalent_radius']) -1 for atom in molecule.atoms])
    a = len(molecule.atoms)
    p = mpc(molecule, 2)
    return float(((a + alpha -1) * ((a + alpha -2)**2)))/((p + alpha)**2)

def s3k(molecule):
    molecule = molecule.hydrogen_suppressed
    alpha = sum([(float(periodic_table[atom.Z]['covalent_radius'])/periodic_table[6]['covalent_radius']) -1 for atom in molecule.atoms])
    a = len(molecule.atoms)
    p = mpc(molecule, 3)
    if a % 2 == 0:
        return float(((a + alpha -3) * ((a + alpha -2)**2)))/((p + alpha)**2)
    else:
        return float(((a + alpha -1) * ((a + alpha -3)**2)))/((p + alpha)**2)

def phi(molecule):
    molecule = molecule.hydrogen_suppressed
    return (s1k(molecule) * s2k(molecule)) / len(molecule.atoms)

def pw(molecule, order):
    molecule = molecule.hydrogen_suppressed
    p = mpc(molecule, order)
    w =  sum(walk_vector(molecule, order).values())/2
    return (float(p)/w)/len(molecule.atoms)

# electrotopologocal index


def intrinsic_state(atom):
    a = (2.0/periodic_table[atom.Z]["principal_quantum_number"]) **2
    b = a * valence_degree(atom) + 1
    c =sum([atom_.Z!= 1 for atom_ in atom.connected_with()])
    return b/c

def intrinsic_state_sum(molecule):
    s = 0
    for atom in molecule.atoms:
        if atom.Z != 1:
            s += intrinsic_state(atom)
    return s

def partition(molecule):
    if cyclomatic_number(molecule) > 0:
        print "It's not acyclic graph"
        return None
    molecule = molecule.hydrogen_suppressed
    p = []
    while len(molecule.atoms) > 1:
        p.append(0)
        if len(molecule.atoms) == 2:
            p[-1] +=2
            molecule.atoms = []
            break
        a = []
        for atom in molecule.atoms:
            if len(atom.bonds) == 1:
                a.append(atom)
        for atom in a:
            atom1 = atom
            bond = atom1.bonds[0]
            for atom_ in bond:
                if atom_ is not atom1:
                    atom2 = atom_
            molecule.atoms.remove(atom1)
            molecule.bonds.remove(bond)
            atom2.bonds.remove(bond)
            p[-1] += 1
    if len(molecule.atoms) > 0:
        if p[-1] <2:
            p[-1] += len(molecule.atoms)
        else:
            p.append(1)
    return p

def bac(molecule):
    p = partition(molecule)
    return sum([x*x for x in p])

def loc(molecule):
    a = len(molecule.atoms)
    p = partition(molecule)
    return -1 * sum([(float(nk)/a)* math.log(float(nk)/a, 2) for nk in p])
