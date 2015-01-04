# -*- coding: utf-8 -*-
__author__ = 'a.shehovtsov'
"""

"""

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


def path_find(molecule):
    def dfs(atom):
        used_atoms.append(atom)

        for bond in atom.bonds:
            next_atom = filter(lambda x: x is not atom,  [atom_ for atom_ in bond])[0]
            #print next_atom
            if next_atom != used_atoms[-2] and next_atom != used_atoms[0]:
                bond_stack.append(bond)
                key = [used_atoms[0], next_atom]
                key.sort()
                key = tuple(key)
                value = tuple(bond_stack[:])
                if key not in path_distance_dct:
                    path_distance_dct[key] = value
                else:
                    if len(path_distance_dct[key]) > len(value):
                        path_distance_dct[key] = value
            if next_atom not in used_atoms:
                dfs(next_atom)
        if bond_stack:
            bond_stack.pop()


    path_distance_dct = {}
    for atom in molecule.atoms:
        used_atoms = [atom]
        bond_stack = []
        dfs(atom)
    return path_distance_dct

def cycle_vertex(molecule):
    """
    Find all vertices which belongs to cycles
    :param molecule:
    :return: set(atom1, atom2, atom3)
    Need tests
    """
    def dfs(atom):
        used_atoms.add(atom)
        atoms = atom.connected_with()
        for next_atom in atoms:
            if next_atom not in used_atoms:
                atom_stack.append(atom)
                dfs(next_atom)
            elif next_atom in used_atoms:
                if atom_stack and next_atom in atom_stack and next_atom is not atom_stack[-1]:
                    i = atom_stack.index(next_atom)
                    cycle = atom_stack[i:] + [atom]
                    for cycle_atom in cycle:
                        cycle_atoms.add(cycle_atom)
        if atom_stack:
            atom_stack.pop()
    used_atoms = set()
    atom_stack = []
    cycle_atoms = set()
    for atom in molecule.atoms:
        if atom not in used_atoms:
            dfs(atom)
    return cycle_atoms
