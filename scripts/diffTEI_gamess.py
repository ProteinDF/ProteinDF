#!/usr/bin/env python

import getopt, sys
import TlGamess
import TlPdfTEI

def main():
    try:
        optlist, args = getopt.gnu_getopt(sys.argv[1:], "hg:p:t:", longopts=["help", "gamess=", "proteindf=", "thresold="])
    except getopt.GetoptError:
        # exit before display the help message
        usage()
        sys.exit(2)


    gamess_out_path_path = ""
    pdf_out_path = ""
    threshold = 1E-9
    for opt, arg in optlist:
        if opt in ("-h", "--help"):
            usage()
            sys.exit(0)
        elif opt in ("-g", "--gaussian"):
            gamess_out_path = arg
        elif opt in ("-p", "--proteindf"):
            pdf_out_path = arg
        elif opt in ("-t", "--threshold"):
            threshold = float(arg)

    #if (len(args) == 0):
        #usage()
        #sys.exit(2)

    print "gamess_out_path =", gamess_out_path
    print "pdf_out_path    =", pdf_out_path
    print "threshold       =", threshold

    # loading gaussian output file
    gamess_parser = TlGamess.GamessParser(gamess_out_path);
    gamess_parser.read()

    pdf = TlPdfTEI.TlPdfTEI(pdf_out_path)
    pdf.read()

    number_of_basis = gamess_parser.get_number_of_basis()
    print "number of basis(GAMESS):", number_of_basis
    
    for i in range(0, number_of_basis):
        for j in range(0, number_of_basis):
            for k in range(0, number_of_basis):
                if (k == i):
                    max_l = j
                else:
                    max_l = k
                for l in range(0, max_l +1):
                    left_v = gamess_parser.get_twoei_from_key(i, j, k, l)
                    right_v = pdf.get_twoei_from_key(i, j, k, l)
                    res_v = left_v - right_v;
                    if (abs(res_v) > threshold):
                        check = "!"
                    else:
                        check = " "
                    print "(%3d %3d | %3d %3d) = % 16.12f  : % 16.12f ==> % 16.12f %s" % (i+1, j+1, k+1, l+1, left_v, right_v, res_v, check)

def usage():
    print "-g gamess.out -p pdf.out -t threshold"


if __name__ == '__main__':
    main()
    
