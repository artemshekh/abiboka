__author__ = 'a.shehovtsov'
"""

"""

def dfs(molecule):
    ex = set()
    def dfs_(atom):
        ex.add(atom)
        for bond in atom.bonds:
            for atom_ in bond:
                if atom_ not in ex:
                    dfs_(atom_)
    for atom in molecule.atoms:
        dfs_(atom)
    print len(ex), len(molecule.atoms)

def check_cycle(molecule):

    def cycle(parent, atom):
        """
        IT's very slow and unefficient check either atom in cycle or not
        for _ in molecule.atoms:
            visited, cycle = [], [False]
            cycle_(_, _)
        :param parent:
        :param atom:
        :return: set of atoms from molecule that situated in cycle
        """
        visited.append(atom)
        next_iteration_atoms = set()
        for bond in atom.bonds:
            for atom_ in bond:
                if atom_ is not parent:
                    if atom_ is visited[0]:
                        cycle_[0] = True
                    if atom_ not in visited:
                        next_iteration_atoms.add(atom_)
        for next_atom in next_iteration_atoms:
                cycle(atom, next_atom)
    s = set()
    for _ in molecule.atoms:
        visited, cycle_ = [], [False]
        cycle(_, _)
        if cycle_[0]:
            s.add(_)
    return s
