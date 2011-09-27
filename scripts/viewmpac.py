#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
display file formatted by MsgPack.
"""

import sys
import optparse
import msgpack

def main():
    # initialize
    parser = optparse.OptionParser(usage="%prog <FILE>",
                                   version="%prog 1.0")

    (opts, args) = parser.parse_args()

    if (len(args) == 0):
        parser.error("input MsgPack file.")
        sys.exit()

    file_path = args[0]

    f = open(file_path, "rb")
    contents = f.read()
    data = msgpack.unpackb(contents)
    f.close()

    print data

if __name__ == '__main__':
    main()
