#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from optparse import OptionParser

import TlVector

def main():
    usage = "usage: %prog [options] vector_file"
    version = "%prog 0.10"
    parser = OptionParser(usage, version=version)

    parser.add_option(
        "-v", "--verbose",
        action = "store_true",
        default = False,
        help = "Be moderately verbose"
        )

    parser.add_option(
        "-g", "--guess_format",
        action = "store_true",
        dest = "guess_format",
        default = False,
        help = "output guess format of ProteinDF"
        )

    (options, args) = parser.parse_args()

    if not args:
        parser.error("requires vector_file")

    file_path = args.pop(0)
    if (options.verbose == True):
        sys.stderr.write("reading file: %s\n" % (file_path))

    v = TlVector.TlVector()

    if (v.is_loadable(file_path) == True):
        v.load_binary(file_path)
        if (options.guess_format == False):
            v.printout()
        else:
            print "TEXT"
            size = v.get_size()
            print size
            print size
            print "0"
            for x in range(size):
                print " % 8.4f" % (v.get(x)),
                if( (x % 10) == 9):
                    print
    else:
        print "unkown file type. stop."
    

if __name__ == '__main__':
    main()
    
