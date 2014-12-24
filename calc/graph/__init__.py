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
