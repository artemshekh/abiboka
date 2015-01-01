# -*- coding: utf-8 -*-
"""

"""

def cached(descriptor):
    name = descriptor.__name__
    def wrapper(molecule):
        if name not in molecule.descriptor_cache:
            molecule.descriptor_cache[name]  = descriptor(molecule)
        return molecule.descriptor_cache[name]
    return wrapper
