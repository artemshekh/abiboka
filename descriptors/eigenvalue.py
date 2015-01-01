"""

Eigenvalue-based descriptors

"""
import math
from calc.matrixes.matrix import AdjacencyMatrix, Matrix, IdentityMatrix
from descriptors.descriptor_utils import p3_matrix, p2_matrix
from utils.periodic_table import periodic_table


def jacobi(matrix):

    def max_(matrix):
        max_ele = 0
        x, y = 0, 1
        for i, row in enumerate(matrix.matrix):
            for j, value in enumerate(row):
                if i != j and abs(value) > max_ele:
                    max_ele, x, y = abs(value), i, j
        return max_ele, x, y

    def rotate(a, p, k, l):
        n = a.rows()
        aDiff = a.matrix[l][l] - a.matrix[k][k]
        if abs(a.matrix[k][l]) < abs(aDiff)*1.0e-36:
            t = a.matrix[k][l]/aDiff
        else:
            phi = aDiff/(2.0*a.matrix[k][l])
            t = 1.0/(abs(phi) + math.sqrt(phi**2 + 1.0))
            if phi < 0.0:
                t = -t
        c = 1.0/math.sqrt(t**2 + 1.0)
        s = t*c
        tau = s/(1.0 + c)
        temp = a.matrix[k][l]
        a.matrix[k][l] = 0.0
        a.matrix[k][k] -= t*temp
        a.matrix[l][l] += t*temp
        for i in range(k):      # Case of i < k
            temp = a.matrix[i][k]
            a.matrix[i][k] = temp - s*(a.matrix[i][l] + tau*temp)
            a.matrix[i][l] += s*(temp - tau*a.matrix[i][l])
        for i in range(k+1, l):  # Case of k < i < l
            temp = a.matrix[k][i]
            a.matrix[k][i] = temp - s*(a.matrix[i][l] + tau*a.matrix[k][i])
            a.matrix[i][l] += s*(temp - tau*a.matrix[i][l])
        for i in range(l+1, n):  # Case of i > l
            temp = a.matrix[k][i]
            a.matrix[k][i] = temp - s*(a.matrix[l][i] + tau*temp)
            a.matrix[l][i] += s*(temp - tau*a.matrix[l][i])
        for i in range(n):      # Update transformation matrix
            temp = p.matrix[i][k]
            p.matrix[i][k] = temp - s*(p.matrix[i][l] + tau*p.matrix[i][k])
            p.matrix[i][l] += s*(temp - tau*p.matrix[i][l])

    n = matrix.rows()
    maxRot = 5*(n**2)       # Set limit on number of rotations
    p = IdentityMatrix.n(n)

    for i in range(maxRot):# Jacobi rotation loop
        aMax, k, l = max_(matrix)
        if aMax < 1.0e-3:
            return matrix
        rotate(matrix, p, k, l)


def ax1(molecule):
    p1 = AdjacencyMatrix.from_molecule(molecule)
    for x in range(len(molecule.atoms)):
        p1.matrix[x] = [math.sqrt(sum(p1.matrix[x])), math.sqrt(periodic_table[molecule.atoms[x].Z]['vdw_radius'])] + p1.matrix[x]
    p1plus = Matrix(p1.matrix)
    matrix = p1plus * p1plus.transpose()
    _ = jacobi(matrix)
    return max([_.matrix[x][x]for x in range(_.rows())])/2


def ax2(molecule):
    p2 = p2_matrix(molecule)
    for x in range(len(molecule.atoms)):
        p2.matrix[x] = [math.sqrt(sum(p2.matrix[x])), math.sqrt(periodic_table[molecule.atoms[x].Z]['vdw_radius'])] + p2.matrix[x]
    p2plus = Matrix(p2.matrix)
    matrix = p2plus * p2plus.transpose()
    _ = jacobi(matrix)
    return max([_.matrix[x][x]for x in range(_.rows())])/2


def ax3(molecule):
    p3 = p3_matrix(molecule)
    for x in range(len(molecule.atoms)):
        p3.matrix[x] = [math.sqrt(sum(p3.matrix[x])), math.sqrt(periodic_table[molecule.atoms[x].Z]['vdw_radius'])] + p3.matrix[x]
    p3plus = Matrix(p3.matrix)
    matrix = p3plus * p3plus.transpose()
    _ = jacobi(matrix)
    return max([_.matrix[x][x]for x in range(_.rows())])/2



eigenvalue_descriptors = {
    'Ax1': ax1,
    'Ax2': ax2,
    'Ax3':ax3
}