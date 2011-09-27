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
    parser = optparse.OptionParser(usage = "%prog [options] PDB_FILE BRD_FILE",
                                   version = "%prog 1.0")
    parser.add_option("-v", "--verbose", dest = "verbose",
                      action="store_false", default = False,
                      help = "print message")
    (opts, args) = parser.parse_args()
        
    if (len(args) < 2):
        parser.print_help()
        sys.exit(1)

    # setting
    pdb_file_path = args[0]
    output_path = args[1]
    verbose = opts.verbose

    # load PDB file
    if (verbose == True):
        print("reading: %s\n" % (pdb_file_path))
    pdb_obj = BrPdb()
    pdb_obj.load(pdb_file_path)
    #print(pdb_obj)
    atom_group = pdb_obj.get_atom_group()
    #print(atom_group)

    # output DfData as MsgPack
    mpac = msgpack.packb(atom_group.get_raw_data())

    # output file
    if (verbose == True):
        print("writing: %s\n" % (output_path))
    savefile = open(output_path, "wb");
    savefile.write(mpac)
    savefile.close()

    # end

if __name__ == '__main__':
    main()

