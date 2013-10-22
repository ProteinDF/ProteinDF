#!/usr/env python

import struct
import copy

class TlVector:
    _header_struct = "<i"
    #_header_struct = "<l" 
    #_header_struct = "<q" # little endian and int_64 support
    _num_of_size = -1
    _vector = []
    
    def __init__(self, size =0):
        self._num_of_size = size
        self._vector = size * [0.0]

    def copy(self):
        answer = copy.deepcopy(self)
        return answer

    def clear(self):
        self._num_of_size = -1
        self._matrix = []

    def resize(self, new_size):
        new_vector = new_size * [0.0]
        for i in range(min(new_size, self._num_of_size)):
            new_vector[i] = self._vector[i]
        self._vector = new_vector

        self._num_of_size = new_size

    def get_size(self):
        return self._num_of_size

    def get(self, index):
        assert ((0 <= index) and (index < self._num_of_size))
        return self._vector[index]
    
    def set(self, index, value):
        assert ((0 <= index) and (index < self._num_of_size))
        self._vector[index] = value

    def is_loadable(self, file_path):
        fin = open(file_path, "rb")
        size_of_header = struct.calcsize(self._header_struct);
        data = fin.read(size_of_header)
        fin.close()
        header = struct.unpack(self._header_struct, data[0: size_of_header])
        size = header[0]
        if (size >= 0):
            return True
        else:
            return False

    def load_binary(self, file_path):
        data = open(file_path, "rb").read()
        # read header
        start = 0
        size_of_header = struct.calcsize(self._header_struct);
        header = struct.unpack(self._header_struct, data[start: start + size_of_header])
        start += size_of_header
        size = header[0]
        assert size >= 0

        # read contents
        self.resize(size)
        size_of_double = struct.calcsize("d");
        for i in range(self._num_of_size):
            value = struct.unpack("d", data[start: start + size_of_double])
            self.set(i, value[0])
            start += size_of_double

    def printout(self):
        for order in range(0, self._num_of_size, 10):
            print "\n",
            for j in range(order, min(order +10, self._num_of_size)):
                print "   %5d th" % (j +1),
            print "\n",
            for j in range(order, min(order +10, self._num_of_size)):
                print "-----------",
            print  "----\n\n",
            for j in range(order, min(order +10, self._num_of_size)):
                print " %10.6lf" % (self.get(j)),
            print "\n",
  
    def __add__(self, other):
        assert isinstance(other, TlVector)
        assert (self.get_size() == other.get_size())

        answer = self.copy()
        answer += other
        return answer

    def __iadd__(self, other):
        assert isinstance(other, TlVector)
        assert (self.get_size() == other.get_size())
            
        for i in range(0, self.get_size()):
            self._vector[i] += other._vector[i]

        return self

    def __sub__(self, other):
        assert isinstance(other, TlVector)

        answer = self.copy()
        answer -= other
        return answer

    def __isub__(self, other):
        assert isinstance(other, TlVector)
        
        for i in range(0, self.get_size()):
            self._vector[i] -= other._vector[i]

        return self



