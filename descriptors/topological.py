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


TODO Sh index 10 eigenvalue for unsymmetric matrix

"""
import operator
import math
from collections import Counter
from descriptors.vertex_degree import intrinsic_state, valence_degree, cluster_coefficient_vertex
from descriptors.vertex_degree import kupchik_vertex_degree, madan_chemical_degree, perturbation_delta_value
from descriptors.vertex_degree import z_delta_number, kier_hall_electronegativity
from descriptors.connectivity import valence_connectivity_index_1
from descriptors.constitutional import number_of_carbon_atoms
from utils.periodic_table import periodic_table
from calc.matrixes.matrix import AdjacencyMatrix, Matrix
from descriptors.descriptor_utils import walk_vector
from descriptors.matrixes import distance_matrix, adjacency_matrix
from descriptors.walk import mpc, path_sequence_matrix
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


@cached
def first_zagreb_index_by_valence_degree(molecule):
    return sum([valence_degree(atom)**2 for atom in molecule.atoms if atom.Z != 1])


@cached
def first_zagreb_index_by_kupchik_degree(molecule):
    return sum([kupchik_vertex_degree(atom)**2 for atom in molecule.atoms if atom.Z != 1])


@cached
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


@cached
def geometric_narumi_index(molecule):
    simple_narumi_index = narumi_simple_index(molecule)
    molecule = molecule.hydrogen_suppressed
    return math.pow(simple_narumi_index, 1.0/molecule.size)


@cached
def total_structure_connectivity_index(molecule):
    return 1.0/math.sqrt(narumi_simple_index(molecule))


@cached
def pogliani_index(molecule):
    descriptor = 0
    for atom in molecule.atoms:
        if atom.Z != 1:
                descriptor += z_delta_number(atom)
    return descriptor


@cached
def ramification_index_1(molecule):
    descriptor = 0
    for bond in molecule.hydrogen_suppressed.bonds:
        atoms = [atom for atom in bond]
        if len(atoms[0].bonds) == 3 and len(atoms[1].bonds) == 3:
            descriptor += 1
    return descriptor


@cached
def ramification_index_2(molecule):
    descriptor = 0
    molecule = molecule.hydrogen_suppressed
    dist_matrix = distance_matrix(molecule)
    for i, row in enumerate(dist_matrix):
        for j, value in enumerate(row[i+1:]):
            if value == 2:
                v1 = len(molecule.atoms[i].bonds)
                v2 = len(molecule.atoms[j+i+1].bonds)
                if v1 == 3 and v2 == 3:
                    descriptor += 1
    return descriptor


@cached
def benzene_like_index(molecule):
    return valence_connectivity_index_1(molecule)/6


@cached
def polarity_wiener_index(molecule):
    descriptor = 0
    dist_matrix = distance_matrix(molecule)
    for row in dist_matrix:
        print row
        for v in row:
            if v == 3:
                descriptor += 1
    return descriptor/2


@cached
def product_of_row_sums(molecule):
    return reduce(operator.mul, [len(atom.bonds) for atom in molecule.hydrogen_suppressed.atoms])


@cached
def log_product_of_row_sums(molecule):
    return math.log(product_of_row_sums(molecule))


@cached
def mean_square_distance_index(molecule):
    descriptor = 0
    dist_matrix = distance_matrix(molecule)
    n = len(dist_matrix)
    for i, row in enumerate(dist_matrix):
        for j, value in enumerate(row[i+1:]):
            descriptor += value**2
    return math.sqrt(float(descriptor)/(n*(n - 1)))


@cached
def superpendentic_index(molecule):
    adj_matrix = adjacency_matrix(molecule)
    dist_matrix = distance_matrix(molecule)
    pendant_vertexes = []
    for i, row in enumerate(adj_matrix):
        if sum(row) == 1:
            pendant_vertexes.append(i)
    descriptor = 0
    for i, row in enumerate(dist_matrix):
        descriptor += reduce(operator.mul, [row[x] for x in pendant_vertexes])
    return math.sqrt(descriptor)


@cached
def petitjean_shape_index(molecule):
    dist_matrix = distance_matrix(molecule)
    diametr, radius = 0, 100000
    for row in dist_matrix:
        _ = max(row)
        if _ > diametr:
            diametr = _
        if _ < radius:
            radius = _
    print diametr, radius
    return (diametr - radius)/float(diametr)


@cached
def eccentricity(molecule):
    dist_matrix = distance_matrix(molecule)
    return sum([max(row) for row in dist_matrix])


@cached
def average_eccentricity(molecule):
    return eccentricity(molecule)/float(molecule.hydrogen_suppressed.size)


@cached
def eccentric(molecule):
    dist_matrix = distance_matrix(molecule)
    av_ecc = average_eccentricity(molecule)
    return sum([abs(max(row) - av_ecc) for row in dist_matrix])/float(molecule.hydrogen_suppressed.size)


@cached
def average_graph_distance_degree(molecule):
    return sum([sum(row) for row in distance_matrix(molecule)])/float(molecule.hydrogen_suppressed.size)


@cached
def mean_distance_degree_deviation(molecule):
    dist_matrix = distance_matrix(molecule)
    return sum([abs(sum(row) - average_graph_distance_degree(molecule)) for row in dist_matrix])\
        / float(molecule.hydrogen_suppressed.size)


@cached
def unipolarity(molecule):
    return min([sum(row) for row in distance_matrix(molecule)])


@cached
def centralization(molecule):
    return sum([sum(row) for row in distance_matrix(molecule)]) \
        - (molecule.hydrogen_suppressed.size * unipolarity(molecule))


@cached
def variance(molecule):
    return max([sum(row) - unipolarity(molecule) for row in distance_matrix(molecule)])


@cached
def radial_centric_information_index(molecule):
    c = Counter()
    for row in distance_matrix(molecule):
        c[max(row)] += 1
    descriptor = 0
    for k, v in c.iteritems():
        ng = float(v)/molecule.hydrogen_suppressed.size
        descriptor += ng * (math.log(ng, 2))
    return -1 * descriptor


@cached
def schultz_topological_index(molecule):
    ad = adjacency_matrix(molecule) + distance_matrix(molecule)
    vertex_degree = Matrix([[x] for x in [len(atom.bonds) for atom in molecule.hydrogen_suppressed.atoms]])
    descriptor = 0
    _ = ad * vertex_degree
    for row in _:
        descriptor += sum(row)
    return descriptor


@cached
def schultz_topological_index_by_valence_degree(molecule):
    ad = adjacency_matrix(molecule) + distance_matrix(molecule)
    vertex_degree = Matrix([[x] for x in [valence_degree(atom) for atom in molecule.atoms if atom.Z != 1]])
    descriptor = 0
    _ = ad * vertex_degree
    for row in _:
        descriptor += sum(row)
    return descriptor


@cached
def gutman_topological_index(molecule):
    vertex_degree = [len(atom.bonds) for atom in molecule.hydrogen_suppressed.atoms]
    descriptor = 0
    for i, row in enumerate(distance_matrix(molecule)):
        for j, value in enumerate(row):
            descriptor += value * vertex_degree[i]*vertex_degree[j]
    return descriptor/2.0


@cached
def gutman_topological_index_by_valence_degree(molecule):
    vertex_degree = [valence_degree(atom) for atom in molecule.atoms if atom.Z != 1]
    descriptor = 0
    for i, row in enumerate(distance_matrix(molecule)):
        for j, value in enumerate(row):
            descriptor += value * vertex_degree[i]*vertex_degree[j]
    return descriptor/2.0


@cached
def xu_index(molecule):
    a = math.sqrt(molecule.hydrogen_suppressed.size)
    vertex_degree = [len(atom.bonds) for atom in molecule.hydrogen_suppressed.atoms]
    m = distance_matrix(molecule)
    distance_degree = [sum(row) for row in m]
    numerator, denumerator = 0, 0
    for i, v in enumerate(vertex_degree):
        numerator += v * distance_degree[i]**2
        denumerator += v * distance_degree[i]
    return a * math.log(float(numerator)/denumerator)


@cached
def mti_index(molecule):
    """
    MTI' index [Muller, Szymanski et al., 1990b; Mihalic, Nikolic et al., 1992]
    S index
    """
    descriptor = 0
    for row in adjacency_matrix(molecule) * distance_matrix(molecule):
        descriptor += sum(row)
    return descriptor


@cached
def sh_index_1(molecule):
    dist_matrix = distance_matrix(molecule)
    vertex_distance_degrees = [sum(row) for row in dist_matrix]
    atoms = molecule.hydrogen_suppressed.atoms
    descriptor = 0
    for bond in molecule.hydrogen_suppressed.bonds:
        atom1, atom2 = [atom for atom in bond]
        atom1_v = valence_degree(atom1)
        atom2_v = valence_degree(atom2)
        i1 = atoms.index(atom1)
        i2 = atoms.index(atom2)
        vdd_1 = vertex_distance_degrees[i1]
        vdd_2 = vertex_distance_degrees[i2]
        descriptor += float(vdd_1 * vdd_2)/(atom1_v * atom2_v)
    return math.log(descriptor)


@cached
def sh_index_2(molecule):
    dist_matrix = distance_matrix(molecule)
    vertex_distance_degrees = [sum(row) for row in dist_matrix]
    atoms = molecule.hydrogen_suppressed.atoms
    descriptor = 0
    for bond in molecule.hydrogen_suppressed.bonds:
        atom1, atom2 = [atom for atom in bond]
        atom1_v = valence_degree(atom1)
        atom2_v = valence_degree(atom2)
        i1 = atoms.index(atom1)
        i2 = atoms.index(atom2)
        vdd_1 = vertex_distance_degrees[i1]
        vdd_2 = vertex_distance_degrees[i2]
        descriptor += float(atom1_v * atom2_v)/(vdd_1 * vdd_2)
    return math.log(descriptor)


@cached
def sh_index_3(molecule):
    dist_matrix = distance_matrix(molecule)
    vertex_distance_degrees = [sum(row) for row in dist_matrix]
    atoms = molecule.hydrogen_suppressed.atoms
    descriptor = 0
    for bond in molecule.hydrogen_suppressed.bonds:
        atom1, atom2 = [atom for atom in bond]
        atom1_v = valence_degree(atom1)
        atom2_v = valence_degree(atom2)
        i1 = atoms.index(atom1)
        i2 = atoms.index(atom2)
        vdd_1 = vertex_distance_degrees[i1]
        vdd_2 = vertex_distance_degrees[i2]
        descriptor += 1/math.sqrt(atom1_v * atom2_v * vdd_1 * vdd_2)
    return math.log(descriptor)


@cached
def sh_index_4(molecule):
    dist_matrix = distance_matrix(molecule)
    vertex_distance_degrees = [sum(row) for row in dist_matrix]
    atoms = molecule.hydrogen_suppressed.atoms
    descriptor = 0
    for bond in molecule.hydrogen_suppressed.bonds:
        atom1, atom2 = [atom for atom in bond]
        atom1_v = valence_degree(atom1)
        atom2_v = valence_degree(atom2)
        i1 = atoms.index(atom1)
        i2 = atoms.index(atom2)
        vdd_1 = vertex_distance_degrees[i1]
        vdd_2 = vertex_distance_degrees[i2]
        descriptor += 1/math.sqrt(float(atom1_v * atom2_v)/(vdd_1 * vdd_2))
    return math.log(descriptor)


@cached
def sh_index_5(molecule):
    dist_matrix = distance_matrix(molecule)
    vertex_distance_degrees = [sum(row) for row in dist_matrix]
    atoms = molecule.hydrogen_suppressed.atoms
    descriptor = 0
    for bond in molecule.hydrogen_suppressed.bonds:
        atom1, atom2 = [atom for atom in bond]
        atom1_v = valence_degree(atom1)
        atom2_v = valence_degree(atom2)
        i1 = atoms.index(atom1)
        i2 = atoms.index(atom2)
        vdd_1 = vertex_distance_degrees[i1]
        vdd_2 = vertex_distance_degrees[i2]
        descriptor += 1/math.sqrt(atom1_v * atom2_v + vdd_1 * vdd_2)
    return math.log(descriptor)


@cached
def sh_index_6(molecule):
    dist_matrix = distance_matrix(molecule)
    vertex_distance_degrees = [sum(row) for row in dist_matrix]
    atoms = molecule.hydrogen_suppressed.atoms
    descriptor = 0
    for bond in molecule.hydrogen_suppressed.bonds:
        atom1, atom2 = [atom for atom in bond]
        atom1_v = valence_degree(atom1)
        atom2_v = valence_degree(atom2)
        i1 = atoms.index(atom1)
        i2 = atoms.index(atom2)
        vdd_1 = vertex_distance_degrees[i1]
        vdd_2 = vertex_distance_degrees[i2]
        descriptor += (atom1_v * atom2_v + vdd_1 * vdd_2)
    return math.log(descriptor)


@cached
def sh_index_7(molecule):
    dist_matrix = distance_matrix(molecule)
    vertex_distance_degrees = [sum(row) for row in dist_matrix]
    atoms = molecule.hydrogen_suppressed.atoms
    descriptor = 0
    for bond in molecule.hydrogen_suppressed.bonds:
        atom1, atom2 = [atom for atom in bond]
        atom1_v = valence_degree(atom1)
        atom2_v = valence_degree(atom2)
        i1 = atoms.index(atom1)
        i2 = atoms.index(atom2)
        vdd_1 = vertex_distance_degrees[i1]
        vdd_2 = vertex_distance_degrees[i2]
        descriptor += (atom1_v * atom2_v + math.log(vdd_1 * vdd_2))
    return math.log(descriptor)


@cached
def sh_index_8(molecule):
    dist_matrix = distance_matrix(molecule)
    vertex_distance_degrees = Matrix([[sum(row)] for row in dist_matrix])
    atoms = molecule.hydrogen_suppressed.atoms
    valence_degrees = Matrix([[valence_degree(atom)] for atom in atoms])
    descriptor = vertex_distance_degrees.transpose() * valence_degrees
    return math.log(descriptor[0][0])


@cached
def sh_index_9(molecule):
    dist_matrix = distance_matrix(molecule)
    vertex_distance_degrees = Matrix([[sum(row)] for row in dist_matrix])
    atoms = molecule.hydrogen_suppressed.atoms
    valence_degrees = Matrix([[valence_degree(atom)] for atom in atoms])
    sd = vertex_distance_degrees * valence_degrees.transpose()
    descriptor = 0
    for row in sd:
        descriptor += sum(row)
    return math.log(descriptor)


@cached
def sh_index_10(molecule):
    dist_matrix = distance_matrix(molecule)
    vertex_distance_degrees = Matrix([[sum(row)] for row in dist_matrix])
    atoms = molecule.hydrogen_suppressed.atoms
    valence_degrees = Matrix([[valence_degree(atom)] for atom in atoms])
    sd = vertex_distance_degrees * valence_degrees.transpose()
    """

    Eigenvalue of matrix
    """
    raise NotImplementedError


@cached
def sh_index(molecule):
    sh1 = sh_index_1(molecule)
    nc = number_of_carbon_atoms(molecule)
    return nc + math.sqrt(nc)*sh1


def eccentric_connectivity_index(molecule):
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
    descriptor = 0
    for atom in molecule.atoms:
        if atom.Z != 1:
            descriptor += intrinsic_state(atom)
    return descriptor


@cached
def intrinsic_state_vector(molecule):
    vector = []
    for atom in molecule.atoms:
        if atom.Z != 1:
            vector.append(intrinsic_state(atom))
    return vector


def electrotopological_state_index_vector(molecule):
    dist_matrix = distance_matrix(molecule)
    descriptor = []
    is_vector = intrinsic_state_vector(molecule)
    for i, atom in enumerate(molecule.hydrogen_suppressed.atoms):
        numerator = map(lambda x: x - is_vector[i], is_vector)
        denominator = map(lambda x: (x+1)**2,  dist_matrix[i])
        descriptor.append(is_vector[i] + sum(map(lambda x, y: float(x)/y, numerator, denominator)))
    return descriptor


def kier_hall_electronegativity_vector(molecule):
    vector = []
    for atom in molecule.atoms:
        if atom.Z != 1:
            vector.append(kier_hall_electronegativity(atom))
    return vector

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
