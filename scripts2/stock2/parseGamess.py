#!/usr/bin/env python

import getopt, sys
import TlGamess

def main():
    try:
        optlist, args = getopt.gnu_getopt(sys.argv[1:], "hg:p:t:", longopts=["help", "gaussian=", "proteindf=", "thresold="])
    except getopt.GetoptError:
        # exit before display the help message
        usage()
        sys.exit(2)

    for opt, arg in optlist:
        if opt in ("-h", "--help"):
            usage()
            sys.exit(0)
        else:
            usage()
            sys.exit(2)

    if (len(args) == 0):
        usage()
        sys.exit(2)

    # gamess_out_path is set by args
    gamess_out_path = args.pop(0)


    # loading gaussian output file
    gamess_parser = TlGamess.GamessParser(gamess_out_path);
    gamess_parser.read()

    number_of_basis = gamess_parser.get_number_of_basis()
    print "number of basis:", number_of_basis
    
    #for i in range(0, number_of_basis):
    #    for j in range(0, number_of_basis):
    #        for k in range(0, number_of_basis):
    #            if (k == i):
    #                max_l = j
    #            else:
    #                max_l = k
    #            for l in range(0, max_l +1):
    #                gau_v = gau.get_twoei_from_key(i, j, k, l)
    #                pdf_v = pdf.get_twoei_from_key(i, j, k, l)
    #                res_v = gau_v - pdf_v;
    #                if (abs(res_v) > threshold):
    #                    check = "*"
    #                else:
    #                    check = " "
    #                print "(%3d %3d | %3d %3d) = % 16.12f  : % 16.12f ==> % 16.12f %s" % (i+1, j+1, k+1, l+1, gau_v, pdf_v, res_v, check)

def usage():
    print "-g gaussian.out -p proteindf.out -t threshold"


if __name__ == '__main__':
    main()
    
