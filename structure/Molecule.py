from calc.graph import Graph

class Molecule():
    def __init__(self):
        self.atoms = set()
        self.bonds = set()

    def add_atom(self, atom):
        self.atoms.add(atom)

    def add_bond(self, bond):
        self.bonds.add(bond)

    def molecular_mass(self):
        return sum([atom.get_relative_mass() for atom in self.atoms])

    def molecular_graph(self):
        return Graph.Graph.from_molecule(self.atoms, self.bonds)

if __name__ == '__main__':
    print Molecule()