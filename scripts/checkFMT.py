#!/usr/bin/env python

import getopt, sys
import re

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

    #
    re_fmt = re.compile("F(\d)\((.+)\)\s*=\s*(\S+)")


def usage():
    print "to implement!"


if __name__ == '__main__':
    main()
    
