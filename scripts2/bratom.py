#!/usr/bin/env python
# -*- coding: utf-8 -*-

import msgpack
import copy

from brposition import *
from brperiodictable import *

class BrAtom(object):
    """
    >>> a = BrAtom()
    >>> a.set_element('Fe')
    >>> a.get_atomic_number()
    26
    """
    __periodicTable = BrPeriodicTable()

    def __init__(self, rhs = None):
        self.data = {}
        self.data['atomic_number'] = 0
        self.data['name'] = "" # CA, CB, etc...
        self.data['coord'] = BrPosition()
        self.data['charge'] = 0.0
        if (isinstance(rhs, BrAtom) == True):
            self.data = copy.deepcopy(rhs.data)
        elif (isinstance(rhs, dict) == True):
            self.set_by_raw_data(rhs)

    def move_to(self, position):
        self.data['coord'].move_to(position)
        return self


    def shift_by(self, direction):
        self.data['coord'] += direction
        return self


    def get_position(self):
        return self.data['coord']


    def set_position(self, position):
        assert(isinstance(position, BrPosition))
        self.data['coord'] = position


    def set_element(self, element):
        value = 0
        if (isinstance(element, int) == True):
            self.data['atomic_number'] = element
        elif (isinstance(element, str) == True):
            self.data['atomic_number'] = self.__periodicTable.get_atomic_number(element)
        else:
            raise BrInputError, "illegal input."


    def get_atomic_number(self):
        return self.data['atomic_number']
    

    def get_symbol(self):
        answer = self.__periodicTable.get_symbol(self.data['atomic_number'])
        return answer


    def set_name(self, name):
        assert(isinstance(name, str))
        self.data['name'] = name


    def get_name(self):
        return self.data['name']


    def get_charge(self):
        return self.data['charge']
    
    
    def set_charge(self, charge):
        assert(isinstance(charge, float))
        self.data['charge'] = float(charge)

    def get_path(self):
        return self.data['path']

    def update_path(self, path):
        self.data['path'] = path

    def set_by_raw_data(self, data):
        self.data = data
        self.data['coord'] = BrPosition(data['coord'])
        return self


    def get_raw_data(self):
        data = copy.deepcopy(self.data)
        data['coord'] = [self.data['coord'].x(),
                         self.data['coord'].y(),
                         self.data['coord'].z()]
        return data

    def get_str(self):
        answer = "name=%s:%s (% e, % e, % e) % f" % (self.get_name(),
                                                     self.get_symbol(),
                                                     self.data['coord'].x(),
                                                     self.data['coord'].y(),
                                                     self.data['coord'].z(),
                                                     self.get_charge())
        return answer

    def __str__(self):
        return self.get_str()
    
    def __eq__(self, rhs):
        answer = False
        if ((isinstance(rhs, BrAtom) == True) and
            (self.get_atomic_number() == rhs.get_atomic_number()) and
            (self.get_position() == rhs.get_position())):
            answer = True
        #if ((isinstance(rhs, BrAtom) == True) and
        #    (self.get_atomic_number() == rhs.get_atomic_number()) and
        #    (abs(self.get_charge() - rhs.get_charge()) < 1.0E-5) and
        #    (self.get_position() == rhs.get_position())):
        #    answer = True
        return answer

    def __ne__(self, rhs):
        return not(self.__eq__(rhs))

    
if __name__ == "__main__":
    import doctest
    doctest.testmod()
