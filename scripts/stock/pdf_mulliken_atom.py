#!/usr/bin/env python
# -*- coding: utf-8 -*-

import getopt, sys
import TlProteinDF

def main():
    # options
    opt_iteration = 0
    opt_verbose = False

    # check arguments
    try:
        optlist, args = getopt.gnu_getopt(sys.argv[1:], "i:vh", longopts=["iteration=", "verbose", "help"])
    except getopt.GetoptError:
        # exit before display the help message
        usage()
        sys.exit(2)

    for opt, arg in optlist:
        if opt in ("-i", "--iteration"):
            opt_iteration = arg
        if opt in ("-v", "--verbose"):
            opt_verbose = True
        if opt in ("-h", "--help"):
            usage()
            sys.exit(0)

    if (len(args) == 0):
        usage()
        sys.exit(2)

    # setting
    fl_out_std_path = args.pop(0)
    if (opt_verbose == True):
        print "loading : %s" % (fl_out_std_path)

    pdf_parser = TlProteinDF.PdfOutput(fl_out_std_path)

    # print
    iteration = opt_iteration
    if (iteration == 0):
        iteration = pdf_parser.get_last_iteration()
    number_of_atoms = pdf_parser.get_number_of_atoms()
    if (opt_verbose == True):
        print "iteration = %d" % (iteration)
        print "number of atoms = %d" % (number_of_atoms)

    for atom_index in range(1, (number_of_atoms +1)):
        atom_symbol, gross_population, mulliken_population = pdf_parser.get_mulliken_atom_population(atom_index)
        print "%5d %2s % 12.6f" % (atom_index, atom_symbol, mulliken_population)
        
def usage():
    print "%s <fl_Out_Std>" % (sys.argv[0])


if __name__ == '__main__':
    main()
    



