from core.exception.exception import MalformedMatrixException


class Matrix(object):
    def __init__(self, matrix=[]):
        self.matrix = matrix

    def rows(self):
        """
        Return number of matrix row
        :return: int
        """
        return len(self.matrix)

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
        return self.rows() == self.cols()

    @classmethod
    def from_file(cls, file):
        matrix = [row.strip('\n').split(' ') for row in file]
        return cls(matrix)

    @classmethod
    def random_matrix(cls, size):
        import random
        matrix = [[random.random() for x in range(size)] for x in range(size)]
        return cls(matrix)

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
            return False
        for row in xrange(self.rows()):
            for col in xrange(self.cols()):
                if not self.matrix[row][col] == self.matrix[col][row]:
                    return False
        return True


if __name__ == '__main__':
    import os
    os.getcwd()
    os.chdir('../../fixtures')
    import time
    start = time.time()
    for x in range(100):
        a = Matrix.random_matrix(5)
        #print a
    print time.time() - start