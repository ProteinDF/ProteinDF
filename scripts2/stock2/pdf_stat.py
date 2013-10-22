#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import optparse
import msgpack

from dfoutput import *

def main():
    # parse options
    #
    parser = optparse.OptionParser(usage="%prog [options] <mpac file>", version="%prog 1.0")
    # set option
    parser.add_option("-v", "--verbose", dest="verbose",
                      help="verbose",
                      action="store_true", default=False)
    (opts, args) = parser.parse_args()

    #if (len(args) == 0):
    #    parser.print_help()
    #    sys.exit(1)

    # variable
    verbose = opts.verbose
    mpac_path = "results.mpac"
    if (len(args) > 0):
        mpac_path = args.pop(0)

    # load
    if (verbose == True):
        sys.stderr.write("loading %s.\n" % (mpac_path))
    f = open(mpac_path, "rb")
    contents = f.read()
    raw_data = msgpack.unpackb(contents)
    f.close()
    
    dfdata = DfData()
    dfdata.set_raw_data(raw_data)
   
    print(dfdata.stat())


if __name__ == '__main__':
    main()
    



