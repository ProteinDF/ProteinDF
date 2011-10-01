#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import copy
import optparse
import msgpack

from brpdb import *

def main():
    usage = """\
%prog [options] from_BridgeFile to_BrideFile
Convert atom \"X\" to corresponding real atom.
"""

    # initialize

    # parse args
    parser = optparse.OptionParser(usage = usage)
    parser.add_option("-v", "--verbose", dest = "verbose",
                      action="store_false", default = False,
                      help = "print message")
    (opts, args) = parser.parse_args()
        
    if (len(args) < 2):
        parser.print_help()
        sys.exit(1)

    # setting
    in_file_path = args[0]
    out_file_path = args[1]
    verbose = opts.verbose

    # reading
    if (verbose == True):
        print("reading: %s\n" % (in_file_path))
    in_file = open(in_file_path, "rb")
    mpac_data = msgpack.unpackb(in_file.read())
    in_file.close()
    
    # prepare atomgroup
    atom_group = BrAtomGroup(mpac_data)
    #print(atom_group)

    atom_group = transX2Ion(atom_group)
    mpac = msgpack.packb(atom_group.get_raw_data())

    # output file
    out_file = open(out_file_path, "wb");
    out_file.write(mpac)
    out_file.close()
    
    # end

def transX2Ion(atom_group):
    assert(isinstance(atom_group, BrAtomGroup) == True)
    new_atom_group = copy.deepcopy(atom_group)

    resname = atom_group.get_name()
    for key in new_atom_group.get_group_list():
        tmp = new_atom_group.get_group(key)
        tmp = transX2Ion(tmp)
        new_atom_group.set_group(key, tmp)

    for key in new_atom_group.get_atom_list():
        atom = new_atom_group.get_atom(key)
        symbol = atom.get_symbol()
        transform_flag = 0
        if (symbol == "X"):
            charge = float(atom.get_charge())
            if (charge > 0.0):
                transform_flag = 1
            elif (charge < 0.0):
                transform_flag = -1
            else:
                if ((resname.upper() == "GLU") or
                    (resname.upper() == "ASP")):
                    transform_flag = 1
                elif ((resname.upper() == "ARG") or
                      (resname.upper() == "LYS")):
                    transform_flag = -1
                else:
                    print("unknown charge: %s @%s" % (atom.get_name(), resname))
        if (transform_flag == 1):
            atom.set_element("Na")
            atom.set_name(" NA ")
            new_atom_group.set_atom(key, atom)
        elif (transform_flag == -1):
            atom.set_element("Cl")
            atom.set_name(" CL ")
            new_atom_group.set_atom(key, atom)
            
    return new_atom_group


if __name__ == '__main__':
    main()

