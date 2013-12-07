#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
convert YAML file to MsgPack.
"""

import sys
import optparse
import msgpack
import yaml

def main():
    usage = """\
%prog [options] YAML_FILE MPAC_FILE
Convert a yaml-formated file to message-pack file.
"""
    # initialize
    parser = optparse.OptionParser(usage = usage)

    (opts, args) = parser.parse_args()

    if (len(args) < 2):
        parser.error("input YAML and output MsgPack file.")
        sys.exit()

    yaml_path = args[0]
    mpac_path = args[1]

    fin = open(yaml_path, "rb")
    contents = fin.read()
    fin.close()

    contents = contents.decode('utf8') # for Japanese
    yaml_data = yaml.load(contents)

    mpac_data = msgpack.packb(yaml_data)

    fout = open(mpac_path, "wb")
    fout.write(mpac_data)
    fout.close()

if __name__ == '__main__':
    main()
