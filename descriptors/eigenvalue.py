"""

Eigenvalue-based descriptors

"""
import math
from calc.matrixes.matrix import AdjacencyMatrix, Matrix, IdentityMatrix
from descriptors.descriptor_utils import p3_matrix, p2_matrix
from utils.periodic_table import periodic_table


def jacobi_rotation(matrix, epsilon= 0.0005):

    def max_element(matrix):
        maximal_value, k, l = 0, None, None
        for i, row in enumerate(matrix.matrix):
            for j, value in enumerate(row[i+1:]):
                if abs(value) > maximal_value and abs(value) > epsilon:
                    maximal_value, k, l = value, i, j+i+1
        return maximal_value, k, l

    def find_phi(matrix, k, l):
        if matrix.matrix[k][k] != matrix.matrix[l][l]:
            angle = float(2*matrix.matrix[k][l])/(matrix.matrix[k][k]-matrix.matrix[l][l])
            phi = 0.5*(math.atan(angle))
        else:
            phi = math.pi/4
        return phi

    def transformation_matrix(n, phi, k, l):
        i = IdentityMatrix.n(n)
        i.matrix[k][l] = -math.sin(phi)
        i.matrix[l][k] = math.sin(phi)
        i.matrix[k][k] = math.cos(phi)
        i.matrix[l][l] = math.cos(phi)
        return i


    n = len(matrix.matrix)
    maxRot = 5*(n**2)
    for i in range(maxRot):
        max, k, l = max_element(matrix)
        if abs(max) < epsilon:
            return matrix
        phi = find_phi(matrix, k, l)
        rt = transformation_matrix(n, phi, k, l)
        tr_transpose = rt.transpose()
        matrix = tr_transpose * (matrix * rt)


def ax1(molecule):
    p1 = AdjacencyMatrix.from_molecule(molecule)
    for x in range(len(molecule.atoms)):
        p1.matrix[x] = [math.sqrt(sum(p1.matrix[x])), math.sqrt(periodic_table[molecule.atoms[x].Z]['vdw_radius'])] + p1.matrix[x]
    p1plus = Matrix(p1.matrix)
    matrix = p1plus * p1plus.transpose()
    _ = jacobi_rotation(matrix, epsilon=0.003)
    return max([_.matrix[x][x]for x in range(_.rows())])/2


def ax2(molecule):
    p2 = p2_matrix(molecule)
    for x in range(len(molecule.atoms)):
        p2.matrix[x] = [math.sqrt(sum(p2.matrix[x])), math.sqrt(periodic_table[molecule.atoms[x].Z]['vdw_radius'])] + p2.matrix[x]
    p2plus = Matrix(p2.matrix)
    matrix = p2plus * p2plus.transpose()
    _ = jacobi_rotation(matrix, epsilon=0.003)
    return max([_.matrix[x][x]for x in range(_.rows())])/2


def ax3(molecule):
    p3 = p3_matrix(molecule)
    for x in range(len(molecule.atoms)):
        p3.matrix[x] = [math.sqrt(sum(p3.matrix[x])), math.sqrt(periodic_table[molecule.atoms[x].Z]['vdw_radius'])] + p3.matrix[x]
    p3plus = Matrix(p3.matrix)
    matrix = p3plus * p3plus.transpose()
    _ = jacobi_rotation(matrix, epsilon=0.003)
    return max([_.matrix[x][x]for x in range(_.rows())])/2



eigenvalue_descriptors = {
    'Ax1': ax1,
    'Ax2': ax2,
    'Ax3':ax3
}