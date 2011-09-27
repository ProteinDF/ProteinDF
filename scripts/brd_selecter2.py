#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import optparse

from bratom import *
from bratomgroup import *
from brselect import *
from brpdb import *

def main():
    # initialize

    # parse args
    parser = optparse.OptionParser(usage = "%prog [options] BRD_FILE",
                                   version = "%prog 1.0")
    parser.add_option("-o", "--output", dest = "output",
                      action="store", default = "",
                      help = "output path")
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

    # reading
    if (verbose == True):
        print("reading: %s" % (mpac_file_path))
    mpac_file = open(mpac_file_path, "rb")
    mpac_data = msgpack.unpackb(mpac_file.read())
    mpac_file.close()
    
    # prepare atomgroup
    atom_group = BrAtomGroup(mpac_data)

    selecter_Na = BrSelect_Atom("Na")
    atom_group_Na = atom_group.select(selecter_Na)
    #print(atom_group_Na)

    selecter_Cl = BrSelect_Atom("Cl")
    atom_group_Cl = atom_group.select(selecter_Cl)

    atom_group_ions = copy.deepcopy(atom_group_Na)
    atom_group_ions |= atom_group_Cl

    #print(atom_group_ions)
    ex = copy.deepcopy(atom_group)
    ex ^= atom_group_ions
    #print(ex)

    pdb_ex = BrPdb()
    pdb_ex.set_by_atomgroup(ex)

    pdb_ion = BrPdb()
    pdb_ion.set_by_atomgroup(atom_group_ions)

    start_res_id = 2001
    start_serial = 20001
    for chain_id in atom_group_Na.get_group_list():
        chain = atom_group_Na.get_group(chain_id)
        for resid in chain.get_group_list():
            res = chain.get_group(resid)
            res.set_name("NA ")
            for atomid in res.get_atom_list():
                res.set_atom(start_serial, res.get_atom(atomid))
                start_serial += 1
                res.erase_atom(atomid)
        for resid in chain.get_group_list():
            chain.set_group(start_res_id, chain.get_group(resid))
            start_res_id += 1
            chain.erase_group(resid)

    start_res_id = 3001
    start_serial = 30001
    for chain_id in atom_group_Cl.get_group_list():
        chain = atom_group_Cl.get_group(chain_id)
        for resid in chain.get_group_list():
            res = chain.get_group(resid)
            res.set_name("CL ")
            for atomid in res.get_atom_list():
                res.set_atom(start_serial, res.get_atom(atomid))
                start_serial += 1
                res.erase_atom(atomid)
        for resid in chain.get_group_list():
            chain.set_group(start_res_id, chain.get_group(resid))
            start_res_id += 1
            chain.erase_group(resid)

    pdb_Na = BrPdb()
    pdb_Na.set_by_atomgroup(atom_group_Na)
    pdb_Cl = BrPdb()
    pdb_Cl.set_by_atomgroup(atom_group_Cl)

    print(pdb_ex)
    #print(pdb_ion)
    print(pdb_Na)
    print(pdb_Cl)
    

if __name__ == '__main__':
    main()

