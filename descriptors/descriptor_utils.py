from calc.matrixes.matrix import Matrix
from collections import defaultdict, Counter

def p3_matrix(molecule):
    matrix = [[0 for y in range(len(molecule.atoms))] for x in range(len(molecule.atoms))]
    d = {}
    for i, _ in enumerate(molecule.atoms):
        d[_] = i
    for i, startatom in enumerate(molecule.atoms):
        walked_atom = set([startatom])
        second_iter_atom = set()
        third_iter_atom = set()
        final_atom = set()
        for bond in startatom.bonds:
            for atom in bond:
                if atom not in walked_atom:
                    second_iter_atom.add(atom)
        for second_atom in second_iter_atom:
            walked_atom.add(second_atom)
            for bond in second_atom.bonds:
                for atom in bond:
                    if atom not in walked_atom:
                        third_iter_atom.add(atom)
        for third_atom in third_iter_atom:
            walked_atom.add(third_atom)
            for bond in third_atom.bonds:
                for atom in bond:
                    if atom not in walked_atom:
                        final_atom.add(atom)
        for atom in final_atom:
            j = d[atom]
            matrix[i][j], matrix[j][i] = 3, 3
    return Matrix(matrix)


def p2_matrix(molecule):
    matrix = [[0 for y in range(len(molecule.atoms))] for x in range(len(molecule.atoms))]
    d = {}
    for i, _ in enumerate(molecule.atoms):
        d[_] = i
    for i, startatom in enumerate(molecule.atoms):
        walked_atom = set([startatom])
        second_iter_atom = set()
        final_atom = set()
        for bond in startatom.bonds:
            for atom in bond:
                if atom not in walked_atom:
                    second_iter_atom.add(atom)
        for second_atom in second_iter_atom:
            walked_atom.add(second_atom)
            for bond in second_atom.bonds:
                for atom in bond:
                    if atom not in walked_atom:
                        final_atom.add(atom)
        for atom in final_atom:
            j = d[atom]
            matrix[i][j], matrix[j][i] = 2, 2
    return Matrix(matrix)



def path_vector(molecule):
    """
    calculate atomic path number for n
    :param molecule:
    :return:
    """
    molecule = molecule.hydrogen_suppressed()
    def dfs(atom):
        used_atom.add(atom)
        bonds = atom.bonds
        for bond in bonds:
            if bond not in bonds_stack:
                for _ in bond:
                    if _ is not atom and _ is not atom_[0]:
                        bonds_stack.append(bond)
                        dct[atom_[0]][len(bonds_stack)] += 1
                        dfs(_)
        if bonds_stack:
            bonds_stack.pop()
    dct = defaultdict(lambda : Counter())
    for atom in molecule.atoms:
        atom_ = [atom]
        used_atom = set()
        bonds_stack = []
        dfs(atom)
    return dct


