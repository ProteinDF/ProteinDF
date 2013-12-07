#!/usr/bin/env python

#import TlPosition

class TlAtom:
    # member variable
    _atom_index = 0
    _atom_symbol = ""
    _position = None
    _charge   = 0.0

    # member function ==================================================
    def __init__(self, atom_index, atom_symbol):
        self._atom_index = atom_index
        self._atom_symbol = atom_symbol

    def get_atom_index(self):
        return self._atom_index

    def get_atom_symbol(self):
        return self._atom_symbol

    def set_charge(self, charge):
        self._charge = charge

    def get_charge(self):
        return self._charge

    
