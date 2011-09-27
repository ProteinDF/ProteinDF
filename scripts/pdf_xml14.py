#!/usr/bin/env python
# -*- coding: utf-8 -*-

import getopt, sys
import xml.dom.minidom
import TlProteinDF14

def main():
    # option
    opt_verbose = False
    
    try:
        optlist, args = getopt.gnu_getopt(sys.argv[1:], "hvs:", longopts=["help", "verbose"])
    except getopt.GetoptError:
        # exit before display the help message
        usage()
        sys.exit(2)

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

    pdf_output = TlProteinDF14.PdfOutput(fl_out_std_path)

    pdf_xml = TlProteinDF14.PdfXml()
    pdf_xml.set_data_from_pdf_output(pdf_output)
    print pdf_xml

def usage():
    print "%s <fl_Out_Std>" % (sys.argv[0])


if __name__ == '__main__':
    main()
    



