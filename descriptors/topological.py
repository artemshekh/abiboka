# -*- coding: utf-8 -*-
"""
Topological descriptors
ZM1 - first Zagreb index (H - depleted)
ZM1_h first Zagreb index (all molecule)
platt_number Platt number
connection_number = connection number (N2)
First Zagreb index by valence vertex degrees
First Zagreb index by Kupchik vertex degrees
First Zagreb index by Madan vertex degrees
Second Zagreb index
Second Zagreb index by valence vertex degrees
Second Zagreb index by Kupchik vertex degrees
Second Zagreb index by Madan vertex degrees

"""
import operator
import math
from collections import Counter
from descriptors.vertex_degree import intrinsic_state, valence_electrones, valence_degree, cluster_coefficient_vertex
from descriptors.vertex_degree import kupchik_vertex_degree, madan_chemical_degree, perturbation_delta_value
from utils.periodic_table import periodic_table
from calc.matrixes.matrix import AdjacencyMatrix, Matrix
from descriptors.descriptor_utils import path_sequence_matrix, walk_vector
from descriptors.walk import mpc
from descriptors.ring_descriptor import cyclomatic_number
from utils.functional import cached


@cached
def first_zagreb_index(molecule):
    """
    First zagreb index
    :param molecule:
    :return: int
    """
    return sum([atom.vertex_degree**2 for atom in molecule.hydrogen_suppressed.atoms])


@cached
def platt_number(molecule):
    """
    Platt number
    :param molecule:
    :return: int
    """
    return first_zagreb_index(molecule) - 2*(molecule.hydrogen_suppressed.size - 1)


@cached
def connection_number(molecule):
    """
    Connection_number
    :param molecule:
    :return: int
    """
    return first_zagreb_index(molecule)/2 - molecule.hydrogen_suppressed.size + 1

close_shell = [1, 3, 11, 19, 37, 55, 87]


def first_zagreb_index_by_valence_degree(molecule):
    return sum([valence_degree(atom)**2 for atom in molecule.atoms if atom.Z != 1])


def first_zagreb_index_by_kupchik_degree(molecule):
    return sum([kupchik_vertex_degree(atom)**2 for atom in molecule.atoms if atom.Z != 1])


def first_zagreb_index_by_madan_degree(molecule):
    return sum([madan_chemical_degree(atom) ** 2 for atom in molecule.atoms if atom.Z != 1])


"""
Need refactor

def permutation_additive(atom, perm_coefficient):
    return valence_degree(atom) + sum([perm_coefficient * valence_degree(atom) for _ in atom.connected_with()])


def permutation_multiplicative(atom, perm_coefficient):
    return valence_degree(atom) + reduce(operator.mul, [perm_coefficient * valence_degree(atom) for _ in atom.connected_with()])

def zm1per(molecule, permutation_coefficient):
    return sum([permutation_additive(atom, permutation_coefficient) ** 2 for atom in molecule.atoms])

def zm1mulper(molecule, permutation_coefficient):
    return sum([permutation_multiplicative(atom, permutation_coefficient) ** 2 for atom in molecule.atoms])
"""


@cached
def second_zagreb_index(molecule):
    _ = []
    for bond in molecule.hydrogen_suppressed.bonds:
        atoms = [atom for atom in bond]
        _.append(atoms[0].vertex_degree * atoms[1].vertex_degree)
    return sum(_)


@cached
def second_zagreb_index_by_valence_degree(molecule):
    _ = []
    for bond in molecule.bonds:
        atoms = [atom for atom in bond]
        _.append(valence_degree(atoms[0]) * valence_degree(atoms[1]))
    return sum(_)


@cached
def second_zagreb_index_by_kupchik_degree(molecule):
    _ = []
    for bond in molecule.bonds:
        atoms = [atom for atom in bond]
        _.append(kupchik_vertex_degree(atoms[0]) * kupchik_vertex_degree(atoms[1]))
    return sum(_)


@cached
def second_zagreb_index_by_madan_degree(molecule):
    _ = []
    for bond in molecule.bonds:
        atoms = [atom for atom in bond]
        _.append(madan_chemical_degree(atoms[0]) * madan_chemical_degree(atoms[1]))
    return sum(_)
"""
Fix
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
"""


@cached
def overall_modified_zagreb_index_0(molecule):
    molecule = molecule.hydrogen_suppressed
    return sum([1.0/len(atom.bonds) for atom in molecule.atoms])


@cached
def overall_modified_zagreb_ondex_by_valence_degree_0(molecule):
    descriptor = []
    for atom in molecule.atoms:
        if atom.Z != 1:
            v = valence_degree(atom)
            if v:
                descriptor.append(1.0/v)
    return sum(descriptor)


@cached
def overall_modified_zagreb_index_1(molecule):
    molecule = molecule.hydrogen_suppressed
    descriptor = []
    for bond in molecule.bonds:
        atoms = [atom for atom in bond]
        descriptor.append(1.0/len(atoms[0].bonds) * 1.0/len(atoms[1].bonds))
    return sum(descriptor)


@cached
def overall_modified_zagreb_index_by_valence_degree_1(molecule):
    descriptor = []
    for bond in molecule.bonds:
        atoms = [atom for atom in bond]
        if all(map(lambda x: x.Z != 1, atoms)):
            descriptor.append(1.0/valence_degree(atoms[0])*1.0/valence_degree(atoms[1]))
    return sum(descriptor)


@cached
def quadratic_index(molecule):
    molecule = molecule.hydrogen_suppressed
    return 3 - 2*molecule.size + first_zagreb_index(molecule)/2.0


@cached
def bertz_branching_index(molecule):
    molecule = molecule.hydrogen_suppressed
    return sum([len(atom.bonds)*(len(atom.bonds) - 1) for atom in molecule.atoms])/2.0


@cached
def narumi_simple_index(molecule):
    molecule = molecule.hydrogen_suppressed
    return math.log(reduce(operator.mul, [len(atom.bonds) for atom in molecule.atoms]))


@cached
def arithmetic_topological_index(molecule):
    molecule = molecule.hydrogen_suppressed
    return 2.0*len(molecule.bonds)/len(molecule.atoms)


@cached
def harmonic_narumi_index(molecule):
    molecule = molecule.hydrogen_suppressed
    return molecule.size/sum([1.0/len(atom.bonds) for atom in molecule.atoms])

def gnar(molecule):
    molecule = molecule.hydrogen_suppressed
    return math.pow(snar(molecule), 1.0/molecule.size)

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
    return math.sqrt(float(_)/(molecule.size*(molecule.size - 1)))

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
    return sum([max(row) for row in m])/float(molecule.size)

def decc(molecule):
    molecule = molecule.hydrogen_suppressed
    m = distance_matrix(molecule)
    aecc_ = sum([max(row) for row in m])/float(molecule.size)
    return sum([abs(max(row) - aecc_) for row in m])/float(molecule.size)

def agdd(molecule):
    molecule = molecule.hydrogen_suppressed
    m = distance_matrix(molecule)
    return sum([sum(row) for row in m])/float(molecule.size)

def mddd(molecule):
    molecule = molecule.hydrogen_suppressed
    m = distance_matrix(molecule)
    agdd_ = sum([sum(row) for row in m])/float(molecule.size)
    return sum([abs(sum(row) - agdd_) for row in m])/float(molecule.size)

def unip(molecule):
    molecule = molecule.hydrogen_suppressed
    m = distance_matrix(molecule)
    return min([sum(row) for row in m])

def cent(molecule):
    molecule = molecule.hydrogen_suppressed
    m = distance_matrix(molecule)
    unip_ = min([sum(row) for row in m])
    return sum([sum(row) for row in m]) - (molecule.size * unip_)

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
        ng = float(v)/molecule.size
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
    a = math.sqrt(molecule.size)
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
    return descriptor/2 + molecule.size

def s1k(molecule):
    molecule = molecule.hydrogen_suppressed
    alpha = sum([(float(periodic_table[atom.Z]['covalent_radius'])/periodic_table[6]['covalent_radius']) -1 for atom in molecule.atoms])
    a = molecule.size
    p = len(molecule.bonds)/2
    return float(((a + alpha) * ((a + alpha -1)**2)))/((p + alpha)**2)

def s2k(molecule):
    molecule = molecule.hydrogen_suppressed
    alpha = sum([(float(periodic_table[atom.Z]['covalent_radius'])/periodic_table[6]['covalent_radius']) -1 for atom in molecule.atoms])
    a = molecule.size
    p = mpc(molecule, 2)
    return float(((a + alpha -1) * ((a + alpha -2)**2)))/((p + alpha)**2)

def s3k(molecule):
    molecule = molecule.hydrogen_suppressed
    alpha = sum([(float(periodic_table[atom.Z]['covalent_radius'])/periodic_table[6]['covalent_radius']) -1 for atom in molecule.atoms])
    a = molecule.size
    p = mpc(molecule, 3)
    if a % 2 == 0:
        return float(((a + alpha -3) * ((a + alpha -2)**2)))/((p + alpha)**2)
    else:
        return float(((a + alpha -1) * ((a + alpha -3)**2)))/((p + alpha)**2)

def phi(molecule):
    molecule = molecule.hydrogen_suppressed
    return (s1k(molecule) * s2k(molecule)) / molecule.size

def pw(molecule, order):
    molecule = molecule.hydrogen_suppressed
    p = mpc(molecule, order)
    w =  sum(walk_vector(molecule, order).values())/2
    return (float(p)/w)/molecule.size


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
    while molecule.size > 1:
        p.append(0)
        if molecule.size == 2:
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
    if molecule.size > 0:
        if p[-1] <2:
            p[-1] += molecule.size
        else:
            p.append(1)
    return p

def bac(molecule):
    p = partition(molecule)
    return sum([x*x for x in p])

def loc(molecule):
    a = molecule.size
    p = partition(molecule)
    return -1 * sum([(float(nk)/a)* math.log(float(nk)/a, 2) for nk in p])

def overall_clustering(molecule):
    return sum([cluster_coefficient_vertex(atom) for atom in molecule.hydrogen_suppressed.atoms])