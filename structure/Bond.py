from structure.Atom import Atom
class Bond(set):
    def __init__(self, order=1, atom1=None, atom2=None, cis_trans=None):
        self.order = order
        self.cis_trans = cis_trans
        self._atoms = set()
        if atom1:
            self.add_atom(atom1)
        if atom2:
            self.add_atom(atom2)

    def add_atom(self, atom):
        if len(self) < 2:
            self.add(atom)
            atom.bonds.append(self)



if __name__ == '__main__':
    bond = Bond(order=1, atom1 = Atom(1), atom2 = Atom(1))
    for atom in bond:
        print atom, atom.bonds