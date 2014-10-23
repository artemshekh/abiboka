# -*- coding: utf-8 -*-
from Vertice import Vertice

class Edge():

    def __init__(self, v1=Vertice(), v2=Vertice()):
        self.vertices = (v1,v2)