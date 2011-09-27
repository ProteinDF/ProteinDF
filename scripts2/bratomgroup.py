#!/usr/bin/env python
# -*- coding: utf-8 -*-

import msgpack
import copy

from bratom import *
from brselect import *

class BrAtomGroup(object):
    """
    >>> a = BrAtomGroup()
    
    """
    def __init__(self, rhs = None):
        self.data = {}
        self.data['atoms'] = {}
        self.data['groups'] = {}
        self.data['name'] = ""
        self.data['charge'] = None
        self.data['path'] = ""
        if (isinstance(rhs, BrAtomGroup) == True):
            self.data = copy.deepcopy(rhs.data)
        elif (isinstance(rhs, dict) == True):
            self.set_by_raw_data(rhs)

    def shift_by(self, direction):
        for key, grp in self.data['groups']:
            grp.shift_by(direction)
        for key, atm in self.data['atoms']:
            atm.shift_by(direction)

    def get_number_of_groups(self):
        return len(self.data['groups'])

    def get_number_of_atoms(self):
        return len(self.data['atoms'])

    def get_number_of_all_atoms(self):
        answer = 0
        for key, grp in self.data['groups']:
            answer += grp.get_number_of_atoms()
        answer += len(self.atoms)
        return answer

    def get_group_list(self):
        return self.data['groups'].keys()

    def has_group(self, key):
        key = str(key)
        return self.data['groups'].has_key(key)

    def set_group(self, key, value):
        key = str(key)
        assert(isinstance(value, BrAtomGroup))
        self.data['groups'][key] = copy.deepcopy(value)

    def get_group(self, key):
        key = str(key)
        return self.data['groups'].get(key, None)

    def erase_group(self, key):
        assert(isinstance(key, str))
        self.data['groups'].pop(key, None)

    def get_atom_list(self):
        return self.data['atoms'].keys()

    def has_atom(self, key):
        key = str(key)
        return self.data['atoms'].has_key(key)

    def get_atom(self, key):
        key = str(key)
        return self.data['atoms'].get(key, None)

    def set_atom(self, key, value):
        key = str(key)
        assert(isinstance(value, BrAtom))
        self.data['atoms'][key] = copy.deepcopy(value)

    def erase_atom(self, key):
        assert(isinstance(key, str))
        self.data['atoms'].pop(key, None)

    def set_name(self, name):
        assert(isinstance(name, str))
        self.data['name'] = name

    def get_name(self):
        return self.data.get('name', "")

    def set_charge(self, value):
        assert(isinstance(value, float))
        self.data['charge'] = value

    def get_charge(self):
        #self.data.setdefault('charge', None)
        #if (self.data['charge'] == None):
        total_charge = 0.0
        for key, group in self.data['groups'].items():
            total_charge += group.get_charge()
        for key, atom in self.data['atoms'].items():
            total_charge += atom.get_charge()
        self.data['charge'] = total_charge
        return self.data['charge']

    def get_path(self):
        return self.data['path']

    def merge(self, rhs):
        assert(isinstance(rhs, BrAtomGroup) == True)
        for key, group in rhs.data['groups'].items():
            self.merge_group(key, group)
        for key, atom in rhs.data['atoms'].items():
            self.set_atom(key, copy.deepcopy(atom))

    def merge_group(self, key, group):
        assert(isinstance(key, str) == True)
        assert(isinstance(group, BrAtomGroup) == True)
        if (self.has_group(key) == True):
            self.data['groups'][key].merge(group)
        else:
            self.set_group(key, group)

    def select(self, selecter):
        assert(isinstance(selecter, BrSelect) == True)
        self.update_path(self.get_path())

        if (selecter.is_match(self) == True):
            return copy.deepcopy(self)
        else:
            answer = BrAtomGroup()
            answer.set_name(self.get_name())
            #answer.set_charge(self.get_charge())
            for key, group in self.data['groups'].items():
                tmp = group.select(selecter)
                if ((tmp.get_number_of_groups() != 0) or
                    (tmp.get_number_of_atoms() != 0)):
                    answer.set_group(key, tmp)
            for key, atom in self.data['atoms'].items():
                if (selecter.is_match(atom) == True):
                    answer.set_atom(key, atom)
            return answer

    def __iand__(self, rhs):
        """
        implement of '&=' operator
        """
        assert(isinstance(rhs, BrAtomGroup) == True)
        self.update_path(self.get_path())
        rhs.update_path(rhs.get_path())

        for key, group in rhs.data['groups'].items():
            if (self.has_group(key) != True):
                self.erase_group(key)
            else:
                self.data['groups'][key] &= rhs.data['groups'][key]
                if ((self.data['groups'][key].get_number_of_groups() == 0) and
                    (self.data['groups'][key].get_number_of_atoms() == 0)):
                    self.erase_group(key)

        for key, atom in rhs.data['atoms'].items():
            if (self.has_atom(key) != True):
                self.erase_atom(key)

        return self

    def __ior__(self, rhs):
        """
        implement of '|=' operator
        """
        assert(isinstance(rhs, BrAtomGroup) == True)
        self.update_path(self.get_path())
        rhs.update_path(rhs.get_path())

        self.merge(rhs)
        return self

    def __ixor__(self, rhs):
        """
        implement of '^=' operator
        """
        assert(isinstance(rhs, BrAtomGroup) == True)
        self.update_path(self.get_path())
        rhs.update_path(rhs.get_path())

        for key, group in rhs.data['groups'].items():
            if (self.has_group(key) == True):
                self.data['groups'][key].__ixor__(group)
                if ((self.data['groups'][key].get_number_of_groups() == 0) and
                    (self.data['groups'][key].get_number_of_atoms() == 0)):
                    self.erase_group(key)
            else:
                self.set_group(key, group)

        for key, atom in rhs.data['atoms'].items():
            if (self.has_atom(key) == True):
                self.erase_atom(key)
            else:
                self.set_atom(key, atom)

        return self

    def update_path(self, path):
        self.data['path'] = path
        for key, group in self.data['groups'].items():
            group.update_path("%s/%s" % (path, key))
        for key, atom in self.data['atoms'].items():
            atom.update_path("%s/%s" % (path, key))
    
    def set_by_raw_data(self, data):
        assert(isinstance(data, dict) == True)
        self.data = {}
        self.data.setdefault('groups', {})
        self.data.setdefault('atoms', {})
        self.data.setdefault('path', "")
        if (data.has_key('groups') == True):
            for key in data['groups'].keys():
                tmp_grp = BrAtomGroup()
                tmp_grp.set_by_raw_data(data['groups'][key])
                #self.data['groups'][key] = copy.deepcopy(tmp_grp)
                self.set_group(key, tmp_grp)
        if (data.has_key('atoms') == True):
            for key in data['atoms'].keys():
                tmp_atom = BrAtom(data['atoms'][key])
                #self.data['atoms'][key] = copy.deepcopy(tmp_atom)
                self.set_atom(key, tmp_atom)
        #self.data['name'] = data.get('name', '')
        self.set_name(data.get('name', ''))
        #self.data['charge'] = data.get('charge', 0.0)
        self.set_charge(data.get('charge', 0.0))
        self.update_path(self.get_path())
        return self

    def get_raw_data(self):
        self.update_path(self.get_path())
        data = {}
        if (len(self.data['groups']) > 0):
            data.setdefault('groups', {})
            for key in self.data['groups'].keys():
                data['groups'][key] = self.data['groups'][key].get_raw_data()
        if (len(self.data['atoms']) > 0):
            data.setdefault('atoms', {})
            for key in self.data['atoms'].keys():
                data['atoms'][key] = self.data['atoms'][key].get_raw_data()
        data['name'] = self.data['name']
        data['charge'] = self.data['charge']
        return data

    
    def get_str(self):
        # for sort
        def cmp_grp(lhs, rhs):
            if (isinstance(lhs, str) == True):
                if (lhs.isdigit() == True):
                    lhs = int(lhs)
                else:
                    lhs = sum([ord(x) for x in lhs])
            if (isinstance(rhs, str) == True):
                if (rhs.isdigit() == True):
                    rhs = int(rhs)
                else:
                    rhs = sum([ord(x) for x in rhs])
            return cmp(lhs, rhs)
        
        answer = ""
        for key in sorted(self.data['groups'].keys(), cmp=cmp_grp):
            answer += self.data['groups'][key].get_str()
        for key in sorted(self.data['atoms'].keys(), cmp=cmp_grp):
            answer += "%s/%s: %s\n" % (self.get_path(), key, self.data['atoms'][key].get_str())
        return answer


    def __str__(self):
        self.update_path(self.get_path())
        return self.get_str()
        

    def __getitem__(self, key):
        assert(isinstance(key, str))
        if (self.data.get('groups', {}).has_key(key) == True):
            return self.data['groups'][key]
        else:
            raise KeyError, key


    def __setitem__(self, key, value):
        assert(isinstance(key, str))
        if (isinstance(value, BrAtomGroup) == True):
            self.data['groups'].__setitem__(key, value)
        elif (isinstance(value, BrAtom) == True):
            self.data['atoms'].__setitem__(key, value)
        else:
            raise ValueError

    
if __name__ == "__main__":
    import doctest
    doctest.testmod()
