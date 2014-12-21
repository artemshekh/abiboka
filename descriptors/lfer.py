"""

lfer descriptors

"""
from utils.periodic_table import periodic_table
from core.exception.exception import TabularDataAbsenceException

mcGowan_volume_dct = {
    1: 8.71,
    6: 16.35,
    7: 14.39,
    8: 12.43,
    9: 10.48,
    15: 24.87,
    16: 22.91,
    17: 20.95,
    35: 26.21,
    53: 34.53

}


def mcGowan_volume(molecule):
    v = 0
    for atom in molecule.atoms:
        if atom.charge:
            print 'Absence of information of McgGwan volume for ions'
            raise TabularDataAbsenceException
        try:
            v += mcGowan_volume_dct[atom.Z]
        except KeyError as err:
            print 'Absense of information mcgowan volume for {}'.format(periodic_table[atom.Z]['name'])
            raise TabularDataAbsenceException
    v -= 6.65 * len(molecule.bonds)
    return v