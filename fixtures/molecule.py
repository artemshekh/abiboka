from structure.Atom import Atom
from structure.Bond import Bond
from structure.Molecule import Molecule

def ethanol():
    m = Molecule()
    m.atoms += [Atom(z=6) for x in range(2)]
    m.atoms += [Atom(z=1) for x in range(6)]
    m.atoms += [Atom(z=8) for x in range(1)]
    atoms = m.atoms
    bonds = [(0,1), (0,2),(0,3),(0,4),(1,8),(1,5),(1,6),(7,8)]
    for i,j in bonds:
        bond = Bond(atom1=atoms[i], atom2=atoms[j])
        m.bonds.append(bond)
        atoms[i].bonds.append(bond)
        atoms[j].bonds.append(bond)
    return m