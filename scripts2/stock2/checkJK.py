#!/usr/bin/env python

import getopt, sys
import TlPdfTEI
#import TlProteinDF
import TlMatrix

def main():
    try:
        optlist, args = getopt.gnu_getopt(sys.argv[1:], "hs:", longopts=["help", "size="])
    except getopt.GetoptError:
        # exit before display the help message
        usage()
        sys.exit(2)

    for opt, arg in optlist:
        if opt in ("-h", "--help"):
            usage()
            sys.exit(0)
        if opt in ("-s", "--size"):
            print "size = ", arg

    if (len(args) == 0):
        usage()
        sys.exit(2)

    # loading gaussian output file
    pdf_out_path = args.pop(0)
    pdf_parser = TlPdfTEI.TlPdfTEI(pdf_out_path);
    pdf_parser.read();

    number_of_basis = pdf_parser.get_number_of_basis()
    print "number_of_basis = ", number_of_basis

    mat_j = pdf_parser.get_J_matrix(number_of_basis)
    print "==== PDF J matrix ===="
    mat_j.printout()
    
    mat_k = pdf_parser.get_K_matrix(number_of_basis)
    print "==== PDF K matrix ===="
    mat_k.printout()

    # PJ
    file_PJ = args.pop(0)
    mat_pj = TlMatrix.TlMatrix_Symmetric()
    mat_pj.load_binary(file_PJ)
    print "==== PDF PJ matrix ===="
    mat_pj.printout()
    
    # PK
    file_PK = args.pop(0)
    mat_pk = TlMatrix.TlMatrix_Symmetric()
    mat_pk.load_binary(file_PK)
    print "==== PDF PK matrix ===="
    mat_pk.printout()

    # print out
    mat_pj -= mat_j
    print "==== PJ - J matrix ===="
    mat_pj.printout()
    
    mat_pk -= mat_k
    print "==== PK - K matrix ===="
    mat_pk.printout()

def usage():
    print "to implement!"


if __name__ == '__main__':
    main()
    
