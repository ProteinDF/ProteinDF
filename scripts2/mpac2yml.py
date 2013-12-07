#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
display file formatted by MsgPack using YAML.
"""

import sys
import optparse
import msgpack
import yaml

def main():
    usage = """\
%prog [options] FILE
Convert a message-pack file to yaml-format and output standard output.
    """
    # initialize
    parser = optparse.OptionParser(usage=usage)

    (opts, args) = parser.parse_args()

    if (len(args) == 0):
        parser.error("input MsgPack file.")
        sys.exit()

    file_path = args[0]

    f = open(file_path, "rb")
    contents = f.read()
    data = msgpack.unpackb(contents)
    f.close()

    yaml_str = yaml.dump(data,
                         encoding='utf8',
                         allow_unicode=True,
                         default_flow_style=False)
    print(yaml_str)

if __name__ == '__main__':
    main()
