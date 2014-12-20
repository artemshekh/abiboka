# -*- coding: utf-8 -*-
"""
Simple descriptor of molecule

"""



_ = {

    #molecular weight
    "MW": lambda x: x.molecular_mass(),

    # average molecular weight
    "AMW": lambda x: x.molecular_mass()/len(x.atoms)
}