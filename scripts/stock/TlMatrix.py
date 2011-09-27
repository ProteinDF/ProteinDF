#!/usr/env python

import struct
import copy

class TlMatrix:
    #_header_struct = "<iii"
    #_header_struct = "<ill" 
    _header_struct = "<iqq" # little endian and int_64 support
    _num_of_rows = -1
    _num_of_cols = -1
    _matrix = []
    
    def __init__(self, row =0, col =0):
        self._num_of_rows = row
        self._num_of_cols = col
        self._matrix = [0.0 for i in range(row * col +1)]

    def copy(self):
        answer = copy.deepcopy(self)
        return answer

    def clear(self):
        self._num_of_rows = -1
        self._num_of_cols = -1
        self._matrix = []

    def resize(self, new_rows, new_cols):
        new_matrix = [0.0 for i in range(new_rows * new_cols +1)]
        for i in range(min((new_rows * new_cols), (self._num_of_rows * self._num_of_cols))):
            new_matrix[i] = self._matrix[i]
        self._matrix = new_matrix
        self._num_of_rows = new_rows
        self._num_of_cols = new_cols

    def get_num_of_rows(self):
        return self._num_of_rows

    def get_num_of_cols(self):
        return self._num_of_cols
        
    def get(self, row, col):
        assert ((0 <= row) and (row < self._num_of_rows))
        assert ((0 <= col) and (col < self._num_of_cols))
        index = row + (self._num_of_rows * col)
        return self._matrix[index]
    
    def set(self, row, col, value):
        assert ((0 <= row) and (row < self._num_of_rows))
        assert ((0 <= col) and (col < self._num_of_cols))
        index = row + (self._num_of_rows * col)
        self._matrix[index] = value

    def is_loadable(self, file_path):
        fin = open(file_path, "rb")
        size_of_header = struct.calcsize(self._header_struct);
        data = fin.read(size_of_header)
        fin.close()
        header = struct.unpack(self._header_struct, data[0: size_of_header])
        matrix_type = header[0]
        row = header[1]  
        col = header[2]
        if (matrix_type == 0):
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
        matrix_type = header[0]
        assert matrix_type == 0
        row = header[1]
        col = header[2]

        # read contents
        self.resize(row, col)
        size_of_double = struct.calcsize("d");
        for r in range(row):
            for c in range(col):
                value = struct.unpack("d", data[start: start + size_of_double])
                self.set(r, c, value[0])
                start += size_of_double

    def printout(self):
        for order in range(0, self._num_of_cols, 10):
            print "       ",
            for j in range(order, min(order +10, self._num_of_cols)):
                print "   %5d th" % (j +1),
            print "\n  ---- ",

            for j in range(order, min(order +10, self._num_of_cols)):
                print "-----------",
            print "\n",
            #print "----\n",

            for i in range(0, self._num_of_rows):
                print " %5d  " % (i +1),
                for j in range(order, min(order +10, self._num_of_cols)):
                    print " %10.6lf" % (self.get(i, j)),
                print "\n",
            print "\n\n",

    def __add__(self, other):
        assert isinstance(other, TlMatrix_Symmetric)
        assert (self.get_num_of_rows() == other.get_num_of_rows())
        assert (self.get_num_of_cols() == other.get_num_of_cols())

        answer = self.copy()
        answer += other
        return answer

    def __iadd__(self, other):
        assert isinstance(other, TlMatrix_Symmetric)
        assert (self.get_num_of_rows() == other.get_num_of_rows())
        assert (self.get_num_of_cols() == other.get_num_of_cols())
            
        for i in range(0, len(self._matrix)):
            self._matrix[i] += other._matrix[i]

        return self

    def __sub__(self, other):
        assert isinstance(other, TlMatrix_Symmetric)
        assert (self.get_num_of_rows() == other.get_num_of_rows())
        assert (self.get_num_of_cols() == other.get_num_of_cols())

        answer = self.copy()
        answer -= other
        return answer

    def __isub__(self, other):
        assert isinstance(other, TlMatrix_Symmetric)
        assert (self.get_num_of_rows() == other.get_num_of_rows())
        assert (self.get_num_of_cols() == other.get_num_of_cols())
            
        for i in range(0, len(self._matrix)):
            self._matrix[i] -= other._matrix[i]

        return self

########################################################################
# 

class TlMatrix_Symmetric:
    #_header_struct = "<ill" 
    _header_struct = "<iqq" # little endian and int_64 support
    _num_of_dims = -1
    _matrix = []
    
    def __init__(self, dim =0):
        self._num_of_dims = dim
        self._matrix = [0.0 for i in range(dim + (dim +1) * dim / 2 +1)]

    def copy(self):
        answer = copy.deepcopy(self)
        return answer
        
    def clear(self):
        self._num_of_dims = -1
        self._matrix = []

    def get_num_of_rows(self):
        return self._num_of_dims

    def get_num_of_cols(self):
        return self._num_of_dims
        
    def resize(self, new_dim):
        new_matrix = [0.0 for i in range(new_dim + (new_dim +1) * new_dim / 2 +1)]
        for i in range(min((new_dim + (new_dim +1) * new_dim / 2), (self._num_of_dims + (self._num_of_dims +1) * self._num_of_dims / 2))):
            new_matrix[i] = self._matrix[i]
        self._matrix = new_matrix
        self._num_of_dims = new_dim

    def get(self, row, col):
        assert ((0 <= row) and (row < self._num_of_dims))
        assert ((0 <= col) and (col < self._num_of_dims))

        if (row < col):
            row, col = col, row
        index = row + (2 * self._num_of_dims - (col +1)) * col / 2;
        return self._matrix[index]

    def set(self, row, col, value):
        assert ((0 <= row) and (row < self._num_of_dims))
        assert ((0 <= col) and (col < self._num_of_dims))
        index = row + (2 * self._num_of_dims - (col +1)) * col / 2;
        self._matrix[index] = value

    def is_loadable(self, file_path):
        fin = open(file_path, "rb")
        size_of_header = struct.calcsize(self._header_struct);
        data = fin.read(size_of_header)
        fin.close()
        header = struct.unpack(self._header_struct, data[0: size_of_header])
        matrix_type = header[0]
        row = header[1]  
        col = header[2]
        if (matrix_type == 2):
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
        matrix_type = header[0]
        assert matrix_type == 2
        row = header[1]  
        col = header[2]
        assert row == col
        dim = row
        self.resize(dim)

        size_of_double = struct.calcsize("d");
        for r in range(dim):
            for c in range(r +1):
                value = struct.unpack("d", data[start: start + size_of_double])
                self.set(r, c, value[0])
                start += size_of_double

        
    def printout(self):
        for order in range(0, self._num_of_dims, 10):
            print "       ",
            for j in range(order, min(order +10, self._num_of_dims)):
                print "   %5d th" % (j +1),
            print "\n  ---- ",

            for j in range(order, min(order +10, self._num_of_dims)):
                print "-----------",
            print "\n",
            #print "----\n",

            for i in range(0, self._num_of_dims):
                print " %5d  " % (i +1),
                for j in range(order, min(order +10, self._num_of_dims)):
                    if (j > i):
                        print "    ----   ",
                    else:
                        print " %10.6lf" % (self.get(i, j)),
                        #print " % 10.e" % (self.get(i, j)),
                print "\n",
            print "\n\n",

    def __add__(self, other):
        assert isinstance(other, TlMatrix_Symmetric)
        assert (self.get_num_of_rows() == other.get_num_of_rows())
        assert (self.get_num_of_cols() == other.get_num_of_cols())

        answer = self.copy()
        answer += other
        return answer

    def __iadd__(self, other):
        assert isinstance(other, TlMatrix_Symmetric)
        assert (self.get_num_of_rows() == other.get_num_of_rows())
        assert (self.get_num_of_cols() == other.get_num_of_cols())
            
        for i in range(0, len(self._matrix)):
            self._matrix[i] += other._matrix[i]

        return self

    def __sub__(self, other):
        assert isinstance(other, TlMatrix_Symmetric)
        assert (self.get_num_of_rows() == other.get_num_of_rows())
        assert (self.get_num_of_cols() == other.get_num_of_cols())

        answer = self.copy()
        answer -= other
        return answer

    def __isub__(self, other):
        assert isinstance(other, TlMatrix_Symmetric)
        assert (self.get_num_of_rows() == other.get_num_of_rows())
        assert (self.get_num_of_cols() == other.get_num_of_cols())
            
        for i in range(0, len(self._matrix)):
            self._matrix[i] -= other._matrix[i]

        return self



