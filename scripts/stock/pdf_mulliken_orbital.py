#!/usr/bin/env python
# -*- coding: utf-8 -*-

import getopt, sys
import TlProteinDF

def main():
    try:
        optlist, args = getopt.gnu_getopt(sys.argv[1:], "vhs:", longopts=["help", "size="])
    except getopt.GetoptError:
        # exit before display the help message
        usage()
        sys.exit(2)

    opt_iteration = 0
    opt_verbose = False

    for opt, arg in optlist:
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
    number_of_orbitals = pdf_parser.get_number_of_orbitals()
    if (opt_verbose == True):
        print "iteration = %d" % (iteration)
        print "number of orbitals = %d" % (number_of_orbitals)

    for gto_index in range(1, (number_of_orbitals +1)):
        atom_index, atom_symbol, shell, gross_population = pdf_parser.get_mulliken_orbital_population(gto_index)
        print "%5d %2s %6d %4s % 12.6f" % (int(atom_index), atom_symbol, gto_index, shell, gross_population)
        
def usage():
    print "%s <fl_Out_Std>" % (sys.argv[0])


if __name__ == '__main__':
    main()
    



