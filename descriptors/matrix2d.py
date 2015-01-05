# -*- coding: utf-8 -*-
"""
2D matrix based descriptors

Balaban-like index from adjacency matrix

"""
import math

from calc.matrixes.matrix import Matrix
from descriptors.eigenvalue import jacobi_rotation
from descriptors.ring_descriptor import cyclomatic_number
from utils.functional import cached


def balaban_like_index_from_adjacency_matrix(molecule):
    descriptor = 0
    molecule = molecule.hydrogen_suppressed
    adjacency_matrix = molecule.adjacency_matrix
    print adjacency_matrix
    for i, row in enumerate(adjacency_matrix):
        for j, value in enumerate(row[i:]):
            if value == 1:
                print i,j, sum(adjacency_matrix[i])*sum(adjacency_matrix[j+i])
                descriptor += 1/math.sqrt(sum(adjacency_matrix[i])*sum(adjacency_matrix[j+i]))
    return descriptor/(cyclomatic_number(molecule) + 1)

@cached
def eigenvalue_from_adjacency_matrix(molecule):
    molecule = molecule.hydrogen_suppressed
    adjacency_matrix = molecule.adjacency_matrix
    m = jacobi_rotation(Matrix(adjacency_matrix), 0.03)
    return [m.matrix[i][i] for i in range(len(m.matrix))]

def spectral_positive_sum_from_adjacency_matrix(molecule):
    eigenvalues = eigenvalue_from_adjacency_matrix(molecule)
    return sum(filter(lambda x: x>0, eigenvalues))

def  normalized_spectral_positive_sum_from_adjacency_matrix(molecule):
    return spectral_positive_sum_from_adjacency_matrix(molecule)/molecule.size

def logarithmic_spectral_positive_sum_from_adjacency_matrix(molecule):
    return math.log(spectral_positive_sum_from_adjacency_matrix(molecule))

def lovasz_pelikan_index(molecule):
    return max(eigenvalue_from_adjacency_matrix(molecule))

def normalized_leading_eigenvalue_from_adjacency_matrix(molecule):
    return lovasz_pelikan_index(molecule)/molecule.size

def spectral_diameter_from_adjacency_matrix(molecule):
    return max(eigenvalue_from_adjacency_matrix(molecule)) - min(eigenvalue_from_adjacency_matrix(molecule))

def spectral_absolute_deviation_from_adjacency_matrix(molecule):
    mean = sum(eigenvalue_from_adjacency_matrix(molecule))/molecule.size
    return sum(abs(value - mean) for value in eigenvalue_from_adjacency_matrix(molecule))

def spectral_mean_absolute_deviation_from_adjacency_matrix(molecule):
    return spectral_absolute_deviation_from_adjacency_matrix(molecule)/molecule.size




