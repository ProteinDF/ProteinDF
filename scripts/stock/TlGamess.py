#!/usr/bin/env python

import re
import TlMatrix

class GamessParser:
    # member variables
    _twoei = {}
    _number_of_basis = 0
    # comiled regulae explations

    # ex)
    #   NUMBER OF CARTESIAN GAUSSIAN BASIS FUNCTIONS =   55
    _re_number_of_basis  = re.compile("\s*NUMBER OF CARTESIAN GAUSSIAN BASIS FUNCTIONS\s+=\s+(\d+)")

    # 2 electron integrals
    # ex)
    #  27   9  14   5  1.0      0.015423783   27   9  14   6  1.0     -0.045896858
    _re_two_electron_integrals = re.compile("\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+\.\d+)\s+([-]?\d+\.\d+)")
    

    # member functions
    def __init__(self, file_path):
        self.file_path = file_path

    def set_file_path(self, file_path):
        self.file_path = file_path

    def get_file_path(self):
        return self.file_path

    def get_number_of_basis(self):
        return self._number_of_basis

    def get_twoei(self):
        return self._twoei

    def get_twoei_from_key(self, i, j, k, l):
        key = (i, j, k, l)
        if self._twoei.has_key(key):
            return self._twoei[key]
        key = (i, j, l, k)
        if self._twoei.has_key(key):
            return self._twoei[key]
        key = (j, i, k, l)
        if self._twoei.has_key(key):
            return self._twoei[key]
        key = (j, i, l, k)
        if self._twoei.has_key(key):
            return self._twoei[key]
        key = (k, l, i, j)
        if self._twoei.has_key(key):
            return self._twoei[key]
        key = (k, l, j, i)
        if self._twoei.has_key(key):
            return self._twoei[key]
        key = (l, k, i, j)
        if self._twoei.has_key(key):
            return self._twoei[key]
        key = (l, k, j, i)
        if self._twoei.has_key(key):
            return self._twoei[key]
        return 0.0
    
    def get_J_matrix(self, max_index):
        #matrix_j = [[0.0 for i in range(0, max_index)] for i in range(0, max_index)]
        matrix_j = TlMatrix.TlMatrix_Symmetric(max_index)
        seq_i = range(0, max_index)
        seq_k = range(0, max_index)
        seq_l = range(0, max_index)
        for i in seq_i:
            seq_j = range(0, i +1)
            for j in seq_j:
                for k in seq_k:
                    for l in seq_l:
                        tmp = matrix_j.get(i , j) + self.get_twoei_from_key(i, j, k, l)
                        matrix_j.set(i, j, tmp)
                        #matrix_j[i][j] += self.get_twoei_from_key(i, j, k, l)
        return matrix_j

    def get_K_matrix(self, max_index):
        #matrix_k = [[0.0 for i in range(max_index +1)] for i in range(max_index +1)]
        matrix_k = TlMatrix.TlMatrix_Symmetric(max_index)
        seq_i = range(0, max_index)
        seq_k = range(0, max_index)
        seq_l = range(0, max_index)
        for i in seq_i:
            seq_j = range(0, i +1)
            for j in seq_j:
                for k in seq_k:
                    for l in seq_l:
                        tmp = matrix_k.get(i, j) + (self.get_twoei_from_key(i, k, j, l) * -0.5)
                        matrix_k.set(i, j, tmp) 
        return matrix_k

    # read and parse file
    def read(self):
        fi = open(self.file_path, "r")
        for line in fi.readlines():
            self._parse_number_of_basis(line)
            self._parse_two_electron_integrals(line)
        fi.close()

    def _parse_number_of_basis(self, line):
        matchObj = self._re_number_of_basis.search(line)
        if (matchObj != None):
            num = int(matchObj.group(1))
            self._number_of_basis = num

    def _parse_two_electron_integrals(self, line):
        line = line.rstrip()
        while (len(line) != 0):
            line = line.lstrip()
            matchObj = self._re_two_electron_integrals.search(line)
            if (matchObj != None):
                i = int(matchObj.group(1))
                j = int(matchObj.group(2))
                k = int(matchObj.group(3))
                l = int(matchObj.group(4))
                t = (i -1, j -1, k -1, l -1)
                v = float(matchObj.group(6))
                self._twoei[t] = v
                #print "GAM(%3d %3d | %3d %3d) = % e" % (i, j, k, l, v)
                #print "before:", line
                line = line[matchObj.end(6):]
                #print "after :", line
            else:
                break

