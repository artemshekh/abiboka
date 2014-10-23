# -*- coding: utf-8 -*-
from collections import Counter


class Vertice():

    def __init__(self, edges=Counter()):
        self.edges = edges