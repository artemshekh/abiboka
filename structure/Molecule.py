from calc.graph import Graph
# TODO! Think about consistency of atom, molecule, and bond class
class Molecule():
    def __init__(self):
        self.atoms = []
        self.bonds = []

    def add_atom(self, atom):
        self.atoms.add(atom)

    def add_atoms(self, atoms):
        for atom in atoms:
            self.add_atom(atom)

    def add_bond(self, bond):
        self.bonds.add(bond)

    def add_bonds(self, bonds):
        for bond in bonds:
            self.add_bond(bond)

    def molecular_mass(self):
        return sum([atom.get_relative_mass() for atom in self.atoms])

    def molecular_graph(self):
        return Graph.Graph.from_molecule(self.atoms, self.bonds)

if __name__ == '__main__':
    print Molecule()