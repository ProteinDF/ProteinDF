#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from optparse import OptionParser

from dfmatrix import *

def main():
    usage = "usage: %prog [options] matrix_file"
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
        parser.error("requires matrix_file")

    file_path = args.pop(0)
    if (options.verbose == True):
        sys.stderr.write("reading file: %s\n" % (file_path))

    matrix_box = DfMatrix()
    matrix_sym = DfSymmetricMatrix()

    if (matrix_box.is_loadable(file_path) == True):
        matrix_box.load_binary(file_path)
        if (options.guess_format == False):
            matrix_box.printout()
        else:
            print "TEXT"
            row = matrix_box.get_num_of_rows()
            col = matrix_box.get_num_of_cols()
            print row
            print col
            for x in range(row):
                for y in range(col):
                    print "  %f" %(matrix_box.get(x, y)),
                print "\n",
            print "\n"
            
    elif (matrix_sym.is_loadable(file_path) == True):
        matrix_sym.load_binary(file_path)
        if (options.guess_format == False):
            matrix_sym.printout()
        else:
            print "TEXT"
            row = matrix_box.get_num_of_rows()
            col = matrix_box.get_num_of_cols()
            print row
            print col
            for x in range(row):
                for y in range(col):
                    print "  %f" %(matrix_box.get(x, y)),
                print "\n",
            print "\n"

    else:
        print "unkown file type. stop."
    

if __name__ == '__main__':
    main()
    
