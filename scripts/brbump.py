#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import copy
import optparse
import msgpack

from bratomgroup import *

class BrBump(object):
    def __init__(self):
        """
        create empty object
        """
        self.__data = []

    def check_atom_group(self, group, within = 2.0):
        output = ""
        atom_list = self.get_atom_list(group)
        for i in range(len(atom_list)):
            pos_i = atom_list[i].get_position()
            for j in range(i +1, len(atom_list)):
                pos_j = atom_list[j].get_position()
                d = pos_i.distanceFrom(pos_j)
                if (d < within):
                   output += "[%s] near [%s] in %f\n" % (atom_list[i].get_path(),
                                                         atom_list[j].get_path(),
                                                         d)
        return output

    def get_atom_list(self, group):
        atom_list = []
        for key in group.get_group_list():
            subgrp = group.get_group(key)
            atom_list = atom_list + self.get_atom_list(subgrp)
        for key in group.get_atom_list():
            atom = group.get_atom(key)
            atom_list.append(copy.deepcopy(atom))
        return atom_list

    
def main():
    # initialize

    # parse args
    parser = optparse.OptionParser(usage = "%prog [options] PDB_FILE",
                                   version = "%prog 1.0")
    parser.add_option("-o", "--output", dest = "output_path",
                      help = "PDB output file", metavar = "FILE")
    parser.add_option("-v", "--verbose", dest = "verbose",
                      action="store_false", default = False,
                      help = "print message")
    (opts, args) = parser.parse_args()
        
    if (len(args) == 0):
        parser.print_help()
        sys.exit(1)

    # setting
    file_path = args[0]
    verbose = opts.verbose

    # reading
    if (verbose == True):
        print("reading: %s\n" % (file_path))
    mpac_file = open(file_path, "rb")
    mpac_data = msgpack.unpackb(mpac_file.read())
    mpac_file.close()
    
    # prepare atomgroup
    atom_group = BrAtomGroup(mpac_data)
    #print(atom_group)

    # 
    bump = BrBump()
    result = bump.check_atom_group(atom_group, 1.0)
    print(result)

    # end

if __name__ == '__main__':
    main()

