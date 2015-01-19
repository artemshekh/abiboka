# -*- coding: utf-8 -*-
"""
Some simple statistics metrics
"""

import itertools
import math


def pearson_coefficient(iterable1, iterable2):
    if len(iterable1) != len(iterable2):
        print 'there must be error'
    mean_it1 = sum(iterable1)/len(iterable1)
    mean_it2 = sum(iterable2)/len(iterable2)
    covariance = 0
    for x, y in itertools.izip(iterable1, iterable2):
        covariance += (x - mean_it1) * (y - mean_it2)

    std_x = math.sqrt(sum([(x - mean_it1)**2 for x in iterable1]))
    std_y = math.sqrt(sum([(y - mean_it2)**2 for y in iterable2]))
    return covariance/(std_x * std_y)