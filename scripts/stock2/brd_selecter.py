#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import optparse

from bratomgroup import *
from bratom import *
from brselect import *

def main():
    # initialize

    # parse args
    parser = optparse.OptionParser(usage = "%prog [options] BRD_FILE",
                                   version = "%prog 1.0")
    parser.add_option("-o", "--output", dest = "output",
                      action="store", default = "",
                      help = "output path")
    parser.add_option("-e", "--regex", dest = "regex",
                      action="store", default = "",
                      help = "regular expression")
    parser.add_option("-n", "--name", dest = "name",
                      action="store", default = "",
                      help = "find by name")
    parser.add_option("-v", "--verbose", dest = "verbose",
                      action="store_true", default = False,
                      help = "print message")
    (opts, args) = parser.parse_args()
        
    if (len(args) == 0):
        parser.print_help()
        sys.exit(1)

    # setting
    mpac_file_path = args[0]
    verbose = opts.verbose
    output_path = opts.output
    if (verbose):
        print("output: %s" % (output_path))
    query_name = opts.name
    if (verbose):
        print("query(name) = [%s]" % (query_name))
    query_regex = opts.regex
    if (verbose):
        print("query(regex) = [%s]" % (query_regex))

    # reading
    if (verbose == True):
        print("reading: %s" % (mpac_file_path))
    mpac_file = open(mpac_file_path, "rb")
    mpac_data = msgpack.unpackb(mpac_file.read())
    mpac_file.close()
    
    # prepare atomgroup
    atom_group = BrAtomGroup(mpac_data)

    if (len(query_name) != 0):
        selecter = BrSelect_Name(query_name)
        atom_group_name = atom_group.select(selecter)
        atom_group = atom_group_name
        #print(atom_group_name)

    if (len(query_regex) != 0):
        selecter = BrSelect_Regex(query_regex)
        atom_group_regex = atom_group.select(selecter)
        atom_group = atom_group_regex
        #print(atom_group_regex)

    # output
    if (len(output_path) != 0):
        if (verbose == True):
            print("writing: %s\n" % (output_path))
        mpac_data = msgpack.packb(atom_group.get_raw_data())
        savefile = open(output_path, "wb");
        savefile.write(mpac_data)
        savefile.close()
    else:
        print(atom_group)

    # end
    

if __name__ == '__main__':
    main()

