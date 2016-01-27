# -*- coding: utf-8 -*-
"""

"""
import math

class Quaternion(object):
    def __init__(self, a, b , c, d):
        self.a = a
        self.b = b
        self.c = c
        self.d = d

    def norm(self):
        return math.sqrt(self.a**2 + self.b**2 + self.c**2 + self.d**2)

    def conjudate(self):
        return Quaternion(self.a, -self.b, -self.c, -self.d)

    def normalize(self, tolerance):

        norm2 = self.a**2 + self.b**2 + self.c**2 + self.d**2
        if abs(norm2 - 1.0) > 0.00001:
            norm = math.sqrt(norm2)
            self.a = self.a/norm
            self.b = self.b/norm
            self.c = self.c/norm
            self.d = self.d/norm

        print self.a, self.b, self.c, self.d

    def mult(self, q):
        w1, x1, y1, z1 = self.a, self.b, self.c, self.d
        w2, x2, y2, z2 = q.a, q.b, q.c, q.d
        w = w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2
        x = w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2
        y = w1 * y2 + y1 * w2 + z1 * x2 - x1 * z2
        z = w1 * z2 + z1 * w2 + x1 * y2 - y1 * x2
        return Quaternion(w, x, y, z)