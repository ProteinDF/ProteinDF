#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import optparse

from bratomgroup import *
from bratom import *
from dfdata import *

def main():
    # initialize

    # parse args
    parser = optparse.OptionParser(usage = "%prog [options] BRD_FILE ProteinDF_results",
                                   version = "%prog 1.0")
    parser.add_option("-v", "--verbose", dest = "verbose",
                      action="store_false", default = False,
                      help = "print this message")
    (opts, args) = parser.parse_args()
        
    if (len(args) == 0):
        parser.print_help()
        sys.exit(1)

    # setting
    brd_path = args[0]
    dfdata_path = args[1]
    verbose = opts.verbose

    # reading
    if (verbose == True):
        print("reading brd: %s\n" % (brd_path))
    brd_fh = open(brd_path, "rb")
    brd_data = msgpack.unpackb(brd_fh.read())
    brd_fh.close()
    atom_group = BrAtomGroup()
    atom_group.set_by_raw_data(brd_data)
    #print(atom_group)

    if (verbose == True):
        print("reading ProteinDF results: %s\n" % (dfdata_path))
    dfdata_fh = open(dfdata_path, "rb")
    dfdata_raw = msgpack.unpackb(dfdata_fh.read())
    dfdata_fh.close()
    dfdata = DfData()
    dfdata.set_raw_data(dfdata_raw)
    #print("atoms=%d" % dfdata.get_number_of_atoms())
    #print("itr=%d" % dfdata.get_number_of_iterations())
    #for i in range(dfdata.get_number_of_atoms()):
    #    print(i, " charge=", dfdata.get_mulliken_atom_population(
    #            dfdata.get_number_of_iterations(), i))
    #print(dfdata)
    
    # apply
    df_atoms = dfdata.get_number_of_atoms()
    iteration = dfdata.get_number_of_iterations()
    residues = atom_group.get_group("_")
    for i in range(df_atoms):
        df_atom = dfdata.get_atom(i)
        df_atom_pos = df_atom.get_position()
        df_atom_pos *= 0.5291772108
        df_atom.move_to(df_atom_pos)
        df_charge = dfdata.get_mulliken_atom_population(iteration, i)
        #print(">>>> %s, charge=%f" % (df_atom.get_str(), df_charge))

        is_found = False
        for resid in residues.get_group_list():
            if (is_found == True):
                break
            residue = residues.get_group(resid)
            for atomid in residue.get_atom_list():
                atom = residue.get_atom(atomid)
                #print("checking...", atom.get_str())
                if (atom == df_atom):
                    atom.set_charge(df_charge)
                    residue.set_atom(atomid, atom)
                    is_found = True
                    #print("FOUND:", atom.get_str())
                    break
        #if (is_found == True):
        #    print("[%4d] FOUND:    " % (i))
        #else:
        #    print("[%4d] NOT FOUND:" % (i))

    #print("-" * 80)
    #print(atom_group)

    for resid in residues.get_group_list():
        residue = residues.get_group(resid)
        charge = residue.get_charge()
        #print("[%s] % f" % (residue.get_str(), charge))
        print("%d, %f" % (int(resid), charge))
        
                
    
    # end

if __name__ == '__main__':
    main()

