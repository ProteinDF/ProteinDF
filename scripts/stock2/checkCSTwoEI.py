#!/usr/bin/env python

import getopt, sys
import TlGaussian
import TlPdfTEI

def main():
    try:
        optlist, args = getopt.gnu_getopt(sys.argv[1:], "hg:p:t:", longopts=["help", "gaussian=", "proteindf=", "thresold="])
    except getopt.GetoptError:
        # exit before display the help message
        usage()
        sys.exit(2)


    gau_path = ""
    pdf_path = ""
    threshold = 0.01
    for opt, arg in optlist:
        if opt in ("-h", "--help"):
            usage()
            sys.exit(0)
        elif opt in ("-g", "--gaussian"):
            gau_path = arg
        elif opt in ("-p", "--proteindf"):
            pdf_path = arg
        elif opt in ("-t", "--threshold"):
            threshold = float(arg)

    #if (len(args) == 0):
        #usage()
        #sys.exit(2)

    print "gaussian = ", gau_path
    print "proteindf = ", pdf_path
    print "threshold = ", threshold

    # loading gaussian output file
    gau = TlGaussian.GaussianParser(gau_path);
    gau.read()

    pdf = TlPdfTEI.TlPdfTEI(pdf_path)
    pdf.read()

    index = gau.get_number_of_basis()
    print "number of basis from gaussian = ", index
    
    for i in range(0, index):
        for j in range(0, index):
            for k in range(0, index):
                if (k == i):
                    max_l = j
                else:
                    max_l = k
                for l in range(0, max_l +1):
                    gau_v = gau.get_twoei_from_key(i, j, k, l)
                    pdf_v = pdf.get_twoei_from_key(i, j, k, l)
                    res_v = gau_v - pdf_v;
                    type_name = pdf.get_twoei_type_name_from_key(i, j, k, l)
                    
                    if (abs(res_v) > threshold):
                        check = "*"
                    else:
                        check = " "
                    print "(%3d %3d | %3d %3d) = % 16.12f  : % 16.12f [%s] ==> % 16.12f %s" % (i+1, j+1, k+1, l+1, gau_v, pdf_v, type_name, res_v, check)

def usage():
    print "-g gaussian.out -p proteindf.out -t threshold"


if __name__ == '__main__':
    main()
    
