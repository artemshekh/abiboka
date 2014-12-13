from calc.graph import Graph

class Molecule():
    def __init__(self):
        self.atoms = set()
        self.bonds = set()
        self.molecular_graph = Graph.Graph.from_molecule(self.atoms, self.bonds)

    def molecular_mass(self):
        return sum([atom.get_relative_mass() for atom in self.atoms])

if __name__ == '__main__':
    print Molecule()