#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import optparse

from brpdb import *
from bratomgroup import *
from bratom import *


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
    mpac_file_path = args[0]
    verbose = opts.verbose

    # reading
    if (verbose == True):
        print("reading: %s\n" % (mpac_file_path))
    mpac_file = open(mpac_file_path, "rb")
    mpac_data = msgpack.unpackb(mpac_file.read())
    mpac_file.close()
    
    # prepare atomgroup
    atom_group = BrAtomGroup(mpac_data)
    #print(atom_group)

    # prepare BrPdb object
    pdb_obj = BrPdb()
    pdb_obj.set_by_atomgroup(atom_group)
    
    # output PDB
    print(pdb_obj)

    # end

if __name__ == '__main__':
    main()

