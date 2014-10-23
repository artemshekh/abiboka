class Matrix(object):
    def __init__(self, matrix=[]):
        self.matrix = matrix

    def rows(self):
        pass

    def cols(self):
        pass

    def is_square(self):
        pass

    @classmethod
    def from_file(cls, file):
        matrix = [row.strip('\n').split(' ') for row in file]
        return cls(matrix)

class AdjacencyMatrix(Matrix):
    def __init__(self, matrix):
        super(AdjacencyMatrix, self).__init__(matrix)
        self.is_adjacent()
        self.n = len(self.matrix[0])

    def is_adjacent(self):
        return True


if __name__ == '__main__':
    import os
    os.getcwd()
    os.chdir('../fixtures')
    AdjacencyMatrix.from_file(open('file_matrix','r'))