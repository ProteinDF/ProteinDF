#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import optparse
import math

class BrPosition(object):
    """
    >>> p = BrPosition([0, 1, 2])
    >>> p.x()
    0.0
    >>> p.y()
    1.0
    >>> p.z()
    2.0

    >>> p == BrPosition([0, 1, 2])
    True

    >>> abs(abs(p) - 2.23606) < 1.0E-5
    True

    >>> p.norm()
    >>> p == BrPosition([0, 1/math.sqrt(5), 2/math.sqrt(5)])
    True
    
    >>> p.move_to([3, 4, 5])
    >>> p == BrPosition([3.0, 4.0, 5.0])
    True

    >>> tmp = p * 2.0
    >>> tmp == BrPosition([6.0, 8.0, 10.0])
    True
    
    >>> tmp = -1.0 * p
    >>> tmp == BrPosition([-6.0, -8.0, -10.0])
    True

    >>> a = BrPosition([1, 2, 3])
    >>> b = BrPosition([2, 3, 4])
    >>> a * b
    20.0

    >>> p = a + b
    >>> p == BrPosition([3, 5, 7])
    True

    >>> a += b
    >>> a == BrPosition([3, 5, 7])
    True

    >>> p = a - b
    >>> p == BrPosition([1, 2, 3])
    True

    >>> a -= b
    >>> a == BrPosition([1, 2, 3])
    True
    
    """
    def __init__(self, rhs = [0.0, 0.0, 0.0]):
        self.epsilon = 1.0E-5
        self.position = [0.0, 0.0, 0.0]
        if ((isinstance(rhs, (list, tuple)) == True) and (len(rhs) == 3)):
            self.position[0] = float(rhs[0])
            self.position[1] = float(rhs[1])
            self.position[2] = float(rhs[2])
        elif (isinstance(rhs, BrPosition) == True):
            self.position = rhs.position
        else:
            raise BrInputError, "illegal input"

    def get_list(self):
        return self.position

    def x(self):
        return self.position[0]


    def y(self):
        return self.position[1]


    def z(self):
        return self.position[2]


    def move_to(self, position):
        tmp = BrPosition(position)
        self.position = tmp.position

    def distanceFrom(self, other):
        d2 = 0.0
        for i in range(3):
            tmp = self.position[i] - other.position[i]
            d2 += tmp * tmp
        return math.sqrt(d2)
    

    def __str__(self):
        answer = "[% e, % e, % e]" % (self.position[0], self.position[1], self.position[2])
        return answer


    def __eq__(self, rhs):
        answer = False
        if (isinstance(rhs, BrPosition) == True):
            if (self.distanceFrom(rhs) < self.epsilon):
                answer = True
        return answer

    def __ne__(self, rhs):
        return not(self.__eq__(rhs))


    def __neg__(self):
        return BrPosition([-x for x in self.position])


    def __abs__(self):
        return math.sqrt(sum([x * x for x in self.position]))


    def norm(self):
        n = self.__abs__()
        self.position = [x / n for x in self.position]


    def __add__(self, rhs):
        if (isinstance(rhs, BrPosition) == True):
            return BrPosition([ x + y for x, y in zip(self.position, rhs.position)])
        else:
            raise BrInputError, "illegal input: BrPosition is required."


    __iadd__ = __add__


    def __sub__(self, rhs):
        return self.__add__(-rhs)


    __isub__ = __sub__


    def __sub__(self, rhs):
        return self.__add__(-rhs)


    def __mul__(self, rhs):
        if (isinstance(rhs, BrPosition) == True):
            return sum([x * y for x, y in zip(self.position, rhs.position)])
        elif (isinstance(rhs, float) == True):
            self.position = [x * rhs for x in self.position]
            return self

    __rmul__ = __mul__


    
if __name__ == "__main__":
    import doctest
    doctest.testmod()
