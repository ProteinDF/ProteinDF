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
    parser = optparse.OptionParser(usage = "%prog [options] CIF_FILE",
                                   version = "%prog 1.0")
    parser.add_option("-o", "--output", dest = "output_path",
                      default = "cif.mpac",
                      help = "output file", metavar = "FILE")
    parser.add_option("-v", "--verbose", dest = "verbose",
                      action="store_false", default = False,
                      help = "print message")
    (opts, args) = parser.parse_args()
        
    if (len(args) == 0):
        parser.print_help()
        sys.exit(1)

    # setting
    cif_file_path = args[0]
    output_path = opts.output_pach
    verbose = opts.verbose

    # load PDB file
    if (verbose == True):
        print("reading: %s\n" % (cif_file_path))
    cif_obj = BrComponentsCif()
    cif_obj.load(cif_file_path)
    atom_group = cif_obj.get_atom_group()
    #print(atom_group)

    # output DfData as MsgPack
    mpac = msgpack.packb(atom_group.get_raw_data())

    # output file
    savefile = open(output_path, "wb");
    savefile.write(mpac)
    savefile.close()

    # end

if __name__ == '__main__':
    main()

