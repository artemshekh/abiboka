import itertools
import operator

from core.exception.exception import MalformedMatrixException


class Matrix(object):

    def __init__(self, matrix=None):
        self.matrix = matrix or []
        self.rows_number = len(self.matrix)
        self.cols_number = self.cols()

    ##################### MAIN OPERATIONS ############################
    def __add__(self, other):
        """
        Addiction of Matrix
        :param other: Matrix
        :return: Matrix
        """
        # We cant add matrices with various dimension
        if not (self.cols_number == other.cols_number and self.rows_number == other.rows_number):
            msg = 'Invalid dimensions. You try to add A {}x{} with B {}x{}'.format(
                str(self.rows_number),
                str(self.cols_number),
                str(other.rows_number),
                str(other.cols_number)
            )
            raise ArithmeticError(msg)
        return Matrix([map(operator.add, *row) for row in itertools.izip(self.matrix, other.matrix)])

    def direct_sum(self, other):
        """
        Direct sum of 2 matrices
        :param other:
        :return:
        """
        matrix = [[0 for col in range(self.cols_number + other.cols_number)] for row in range(self.rows_number + other.rows_number)]
        for i, row in enumerate(self.matrix):
            for j, col in enumerate(row):
                matrix[i][j] = col
        for i, row in enumerate(other.matrix):
            for j, col in enumerate(row):
                matrix[i + self.rows_number][j + self.cols_number] = col
        return Matrix(matrix)

    def kroneker_product(self, other):
        """
        Kroneker product of 2 matrix
        :param other:
        :return:
        """
        matrix = [[0 for x in range(self.cols_number * other.cols_number)] for x in range(self.rows_number * other.rows_number)]
        for i1, row1 in enumerate(self.matrix):
            for j1, column1 in enumerate(row1):
                for i2, row2 in enumerate(other.matrix):
                    for j2, column2 in enumerate(row2):
                        matrix[other.rows_number*i1 + i2][other.cols_number*j1 + j2] = column1 * column2
        return Matrix(matrix)

    def kroneker_sum(self, other):
        matrix = self.kroneker_product(Matrix.create_identity(other.rows_number))\
                 + other.kroneker_product(Matrix.create_identity(self.rows_number))
        return matrix

    def __mul__(self, other):
        if not self.cols_number == other.rows_number:
            raise ArithmeticError('number cols of first matrix = number rows of second matrix')
        matrix = [[0 for col in range(other.cols_number)] for row in range(self.rows_number)]
        for rowindex, row in enumerate(matrix):
            for colindex, col in enumerate(row):
                for l, x in enumerate(self.matrix[rowindex]):
                    matrix[rowindex][colindex] += x * other.matrix[l][colindex]
        return Matrix(matrix)

    def cols(self):
        """
        Return number of columns
        :return: int
        """
        var = [len(row) for row in self.matrix]
        if min(var) == max(var):
            return min(var)
        else:
            raise MalformedMatrixException

    def is_square(self):
        """
        Matrix is square
        :return: bool
        """
        return self.rows_number == self.cols_number

    def is_diagonal(self):
        if not self.is_square():
            return False
        for i, row in enumerate(self.matrix):
            for j, col in enumerate(row):
                if i != j and self.matrix[i][j] != 0:
                    return False
        return True

    def is_identity(self):
        if not self.is_square():
            return False
        for i, row in enumerate(self.matrix):
            for j, column in enumerate(row):
                if i != j:
                    if self.matrix[i][j] != 0 or self.matrix[i][j] != self.matrix[j][i]:
                        return False
                else:
                    if self.matrix[i][j] != 1:
                        return False
        return True

    def __eq__(self, other):
        if not (self.cols_number == other.cols_number and self.rows_number == other.rows_number):
            return False
        for row in range(self.rows_number):
            for col in range(self.cols_number):
                if not self.matrix[row][col] == other.matrix[row][col]:
                    return False
        return True


    @classmethod
    def from_file(cls, filename):
        matrix = [row.strip('\n').split(' ') for row in filename]
        return cls(matrix)

    @classmethod
    def random_matrix(cls, size):
        import random
        matrix = [[random.random() for x in range(size)] for x in range(size)]
        return cls(matrix)

    @classmethod
    def create_identity(cls, size):
        matrix = [[0 for x in range(size)] for x in range(size)]
        for i, row in enumerate(matrix):
            for j, column in enumerate(row):
                if i == j:
                    matrix[i][j] = 1
        return Matrix(matrix)

    def __str__(self):
        """
        Matrix string
        :return: str
        """
        return '\n'.join(['|' + ' '.join(map(str, row)) + '|' for row in self.matrix])


class AdjacencyMatrix(Matrix):
    def __init__(self, matrix):
        super(AdjacencyMatrix, self).__init__(matrix)
        self.is_adjacent()
        self.n = len(self.matrix[0])

    def is_adjacent(self):
        if not self.is_square():
            raise MalformedMatrixException('Adjacency matrix must be square')
        for row in xrange(self.rows_number):
            for col in xrange(row, self.cols_number):
                if not self.matrix[row][col] == self.matrix[col][row]:
                    return False
        return True


if __name__ == '__main__':
    import os
    os.getcwd()
    os.chdir('../../fixtures')
    import time
    s = time.time()
    for x in range(800):
        a = Matrix.random_matrix(100)
        b = Matrix.random_matrix(100)
        c = a + b
    print(time.time() - s)