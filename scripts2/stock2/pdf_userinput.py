#!/usr/bin/env python
# -*- coding: utf-8 -*-

import getopt, sys
import PdfUserinput

def main():
    opt_verbose = False

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

    # setting
    fl_globalinput_path = args.pop(0)
    print "loading : %s" % (fl_globalinput_path)

    pdf_userinput = PdfUserinput.PdfUserinput(fl_globalinput_path)

    print pdf_userinput

def usage():
    print "%s <fl_Out_Std>" % (sys.argv[0])


if __name__ == '__main__':
    main()
    



