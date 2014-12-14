from calc.graph import Graph

class Molecule():
    def __init__(self):
        self.atoms = set()
        self.bonds = set()

    def molecular_mass(self):
        return sum([atom.get_relative_mass() for atom in self.atoms])

    def molecular_graph(self):
        return Graph.Graph.from_molecule(self.atoms, self.bonds)

if __name__ == '__main__':
    print Molecule()