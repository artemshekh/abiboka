from structure.Atom import Atom
class Bond():
    def __init__(self, order=1, atom1=None, atom2=None, cis_trans=None):
        self.order = order
        self.cis_trans = cis_trans
        self.atoms = set([atom1, atom2])


if __name__ == '__main__':
    print Bond(order=1, atom1 = Atom(1), atom2 = Atom(1))