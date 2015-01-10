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

    def __hash__(self):
        return id(self)

    def is_aromatic(self):
        return all([atom_.aromatic for atom_ in self])

    def is_amide_bond(self):
        atoms = [atom for atom in self]
        if atoms[0].Z == 6 and atoms[1].Z ==7:
            for new_bond in atoms[0].bonds:
                if new_bond.order ==2 and any(map(lambda x: x.Z == 8, [atom for atom in new_bond])):
                    return True
        if atoms[0].Z == 7 and atoms[1].Z ==6:
            for new_bond in atoms[1].bonds:
                if new_bond.order ==2 and any(map(lambda x: x.Z == 8, [atom for atom in new_bond])):
                    return True
        return False

    def is_terminal(self):
        return any(map(lambda x: x.vertex_degree == 1, [atom for atom in self]))

    def adjacent_to_triple_bond(self):
        atoms = [atom for atom in self]
        for atom in atoms:
            for bond in atom.bonds:
                if bond.order == 3:
                    return True

    @property
    def conventional_bond_order(self):
        if self.is_aromatic():
            return 1.5
        else:
            return self.order


if __name__ == '__main__':
    bond = Bond(order=1, atom1 = Atom(1), atom2 = Atom(1))
    for atom in bond:
        print atom, atom.bonds