#!/usr/env python

import os
import struct
import copy
import array

class DfMatrix(object):
    __header_struct = "<iii"
    #_header_struct = "<ill" 
    #_header_struct = "<iqq" # little endian and int_64 support
    
    def __init__(self, row =0, col =0, raw_data =None):
        self.__num_of_rows = row
        self.__num_of_cols = col
        need_init = True

        # setup from raw_data
        if (isinstance(raw_data, dict) == True):
            self.__num_of_rows = raw_data.get('row', 0)
            self.__num_of_cols = raw_data.get('col', 0)
            matrix_type = raw_data.get('type', None)
            buf = raw_data.get('data', None)
            if ((matrix_type == 'RSFD') and (buf != None)):
                self.__matrix = array.array('d', buf)
                need_init = False

        if (need_init == True):
            self.__matrix = array.array('d', [0.0 for i in xrange(row * col +1)])
                

    def copy(self):
        answer = copy.deepcopy(self)
        return answer

    def clear(self):
        self.__num_of_rows = -1
        self.__num_of_cols = -1
        self.__matrix = array.array('d')

    def resize(self, new_rows, new_cols):
        new_matrix = array.array('d', [0.0 for i in xrange(new_rows * new_cols +1)])
        for i in xrange(min((new_rows * new_cols), (self.get_num_of_rows() * self.get_num_of_cols()))):
            new_matrix[i] = self.__matrix[i]
        self.__matrix = new_matrix
        self.__num_of_rows = new_rows
        self.__num_of_cols = new_cols

    def get_num_of_rows(self):
        return self.__num_of_rows

    def get_num_of_cols(self):
        return self.__num_of_cols
        
    def get(self, row, col):
        assert ((0 <= row) and (row < self.get_num_of_rows()))
        assert ((0 <= col) and (col < self.get_num_of_cols()))
        index = row + (self.get_num_of_rows() * col)
        return self.__matrix[index]
    
    def set(self, row, col, value):
        assert ((0 <= row) and (row < self.get_num_of_rows()))
        assert ((0 <= col) and (col < self.get_num_of_cols()))
        index = row + (self.get_num_of_rows() * col)
        self.__matrix[index] = value

    def is_loadable(self, file_path):
        fin = open(file_path, "rb")
        size_of_header = struct.calcsize(self.__header_struct);
        data = fin.read(size_of_header)
        fin.close()
        header = struct.unpack(self.__header_struct, data[0: size_of_header])
        matrix_type = header[0]
        row = header[1]  
        col = header[2]
        if (matrix_type == 0):
            return True
        else:
            return False

    def load_binary(self, file_path):
        if (os.path.exists(file_path) == True):
            fin = open(file_path, "rb")
            # read header
            size_of_header = struct.calcsize(self.__header_struct);
            header_bin = fin.read(size_of_header)
            header = struct.unpack(self.__header_struct, header_bin)
            matrix_type = header[0]
            row = header[1]
            col = header[2]
            assert(matrix_type == 0)

            # read contents
            self.resize(row, col)
            size_of_double = struct.calcsize("d");
            for r in xrange(row):
                for c in xrange(col):
                    value = struct.unpack("d", fin.read(size_of_double))
                    self.set(r, c, value[0])
            fin.close()
        else:
            print("file not found: %s" % (file_path))

    def printout(self):
        for order in xrange(0, self.get_num_of_cols(), 10):
            print "       ",
            for j in xrange(order, min(order +10, self.get_num_of_cols())):
                print "   %5d th" % (j +1),
            print "\n  ---- ",

            for j in xrange(order, min(order +10, self.get_num_of_cols())):
                print "-----------",
            print "\n",
            #print "----\n",

            for i in xrange(0, self.get_num_of_rows()):
                print " %5d  " % (i +1),
                for j in xrange(order, min(order +10, self.get_num_of_cols())):
                    print " %10.6lf" % (self.get(i, j)),
                print "\n",
            print "\n\n",


    def get_raw_data(self):
        raw = {}
        raw['row'] = self.get_num_of_rows()
        raw['col'] = self.get_num_of_cols()
        raw['type'] = 'RSFD'
        raw['data'] = self.__matrix
        return raw

    def __add__(self, other):
        assert isinstance(other, DfMatrix)
        assert (self.get_num_of_rows() == other.get_num_of_rows())
        assert (self.get_num_of_cols() == other.get_num_of_cols())

        answer = self.copy()
        answer += other
        return answer

    def __iadd__(self, other):
        assert isinstance(other, DfMatrix)
        assert (self.get_num_of_rows() == other.get_num_of_rows())
        assert (self.get_num_of_cols() == other.get_num_of_cols())
            
        for i in xrange(0, len(self.__matrix)):
            self.__matrix[i] += other.__matrix[i]

        return self

    def __sub__(self, other):
        assert isinstance(other, DfMatrix)
        assert (self.get_num_of_rows() == other.get_num_of_rows())
        assert (self.get_num_of_cols() == other.get_num_of_cols())

        answer = self.copy()
        answer -= other
        return answer

    def __isub__(self, other):
        assert isinstance(other, DfMatrix)
        assert (self.get_num_of_rows() == other.get_num_of_rows())
        assert (self.get_num_of_cols() == other.get_num_of_cols())
            
        for i in xrange(0, len(self.__matrix)):
            self.__matrix[i] -= other.matrix[i]

        return self

########################################################################
# 

class DfSymmetricMatrix(object):
    _header_struct = "<ill" 
    #_header_struct = "<iqq" # little endian and int_64 support
    
    def __init__(self, dim =0, raw_data =None):
        self.__num_of_dims = dim
        need_init = True

        # setup from raw_data
        if (isinstance(raw_data, dict) == True):
            self.__num_of_dims = raw_data.get('row', 0)
            assert(self.__num_of_dims == raw_data.get('col', 0))
            matrix_type = raw_data.get('type', None)
            buf = raw_data.get('data', None)
            if ((matrix_type == 'RLHD') and (buf != None)):
                self.__matrix = array.array('d', buf)
                need_init = False

        if (need_init == True):
            self.__matrix = array.array('d', [0.0 for i in xrange(dim + (dim +1) * dim / 2 +1)])

    def copy(self):
        answer = copy.deepcopy(self)
        return answer
        
    def clear(self):
        self.___num_of_dims = -1
        self.__matrix = array.array('d')

    def get_num_of_rows(self):
        return self.__num_of_dims

    def get_num_of_cols(self):
        return self.__num_of_dims
        
    def resize(self, new_dim):
        dim = self.get_num_of_rows()
        assert(dim == self.get_num_of_cols())
        new_matrix = array.array('d', [0.0 for i in xrange(new_dim + (new_dim +1) * new_dim / 2 +1)])
        for i in xrange(min((new_dim + (new_dim +1) * new_dim / 2), (dim + (dim +1) * dim / 2))):
            new_matrix[i] = self.__matrix[i]
        self.__matrix = new_matrix
        self.__num_of_dims = new_dim


    def get(self, row, col):
        assert ((0 <= row) and (row < self.get_num_of_rows()))
        assert ((0 <= col) and (col < self.get_num_of_cols()))
        dim = self.get_num_of_rows()
        assert(dim == self.get_num_of_cols())

        if (row < col):
            row, col = col, row
        index = row + (2 * dim - (col +1)) * col / 2;
        return self.__matrix[index]


    def set(self, row, col, value):
        assert ((0 <= row) and (row < self.get_num_of_rows()))
        assert ((0 <= col) and (col < self.get_num_of_cols()))
        dim = self.get_num_of_rows()
        assert(dim == self.get_num_of_cols())
        index = row + (2 * dim - (col +1)) * col / 2;
        self.__matrix[index] = value


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
        if (os.path.exists(file_path) == True):
            fin = open(file_path, "rb")
            # read header
            size_of_header = struct.calcsize(self._header_struct);
            header = struct.unpack(self._header_struct,
                               fin.read(size_of_header))
            matrix_type = header[0]
            row = header[1]  
            col = header[2]
            dim = row
            assert(matrix_type == 2)
            assert(row == col)
            self.resize(dim)

            size_of_double = struct.calcsize("d");
            for r in xrange(dim):
                for c in xrange(r +1):
                    value = struct.unpack("d", fin.read(size_of_double))
                    self.set(r, c, value[0])
            fin.close()
        else:
            print("file not found: %s" % (file_path))
        
    def printout(self):
        dim = self.get_num_of_rows()
        assert(dim == self.get_num_of_cols())
        for order in xrange(0, dim, 10):
            print "       ",
            for j in xrange(order, min(order +10, dim)):
                print "   %5d th" % (j +1),
            print "\n  ---- ",

            for j in xrange(order, min(order +10, dim)):
                print "-----------",
            print "\n",
            #print "----\n",

            for i in xrange(0, dim):
                print " %5d  " % (i +1),
                for j in xrange(order, min(order +10, dim)):
                    if (j > i):
                        print "    ----   ",
                    else:
                        print " %10.6lf" % (self.get(i, j)),
                        #print " % 10.e" % (self.get(i, j)),
                print "\n",
            print "\n\n",


    def get_raw_data(self):
        raw = {}
        raw['row'] = self.get_num_of_rows()
        raw['col'] = self.get_num_of_cols()
        raw['type'] = 'RSFD'
        raw['data'] = self.__matrix
        return raw


    def __add__(self, other):
        assert isinstance(other, DfSymmetricMatrix)
        assert (self.get_num_of_rows() == other.get_num_of_rows())
        assert (self.get_num_of_cols() == other.get_num_of_cols())

        answer = self.copy()
        answer += other
        return answer

    def __iadd__(self, other):
        assert isinstance(other, DfSymmetricMatrix)
        assert (self.get_num_of_rows() == other.get_num_of_rows())
        assert (self.get_num_of_cols() == other.get_num_of_cols())
            
        for i in xrange(0, len(self.__matrix)):
            self.__matrix[i] += other.__matrix[i]

        return self

    def __sub__(self, other):
        assert isinstance(other, DfSymmetricMatrix)
        assert (self.get_num_of_rows() == other.get_num_of_rows())
        assert (self.get_num_of_cols() == other.get_num_of_cols())

        answer = self.copy()
        answer -= other
        return answer

    def __isub__(self, other):
        assert isinstance(other, DfSymmetricMatrix)
        assert (self.get_num_of_rows() == other.get_num_of_rows())
        assert (self.get_num_of_cols() == other.get_num_of_cols())
            
        for i in xrange(0, len(self.__matrix)):
            self.__matrix[i] -= other.__matrix[i]

        return self



