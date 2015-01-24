# -*- coding: utf-8 -*-
"""

"""


def cached(descriptor):
    name = descriptor.__name__
    def wrapper(molecule, *args):
        if name not in molecule.descriptor_cache:
            molecule.descriptor_cache[name] = descriptor(molecule, *args)
        return molecule.descriptor_cache[name]
    return wrapper
