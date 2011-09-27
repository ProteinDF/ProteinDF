#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
convert the file written by MsgPack to AVS Field data.
"""

import sys
import optparse
import msgpack
import struct

def main():
    usage = """\
%prog [options] <mpac_FILE> <field_data_FILE>
Convert mpac file to AVS field data.
    """
    # initialize
    parser = optparse.OptionParser(usage=usage)

    (opts, args) = parser.parse_args()

    if (len(args) == 0):
        parser.error("input MsgPack file.")
        sys.exit()

    mpac_file_path = args[0]
    fld_file_path = args[1]

    fin = open(mpac_file_path, "rb")
    contents = fin.read()
    data = msgpack.unpackb(contents)
    fin.close()

    # read coord
    num_of_grid_x = data["num_of_grid_x"]
    num_of_grid_y = data["num_of_grid_y"]
    num_of_grid_z = data["num_of_grid_z"]
    coord = data["coord"]

    # read data
    pa = data["MO_34"]

    # prepare
    num_of_grids = num_of_grid_x * num_of_grid_y * num_of_grid_z
    print("num of grids = %d" % (num_of_grids))

    # output
    nspace = 3
    veclen = 1
    data_type = "float"
    field_type = "irregular"
    label = "none"

    fout = open(fld_file_path, "wb")
    fout.write("# AVS field file\n")
    fout.write("ndim = 3\n")
    fout.write("dim1 = %d\n" % (num_of_grid_x))
    fout.write("dim2 = %d\n" % (num_of_grid_y))
    fout.write("dim3 = %d\n" % (num_of_grid_z))
    fout.write("nspace = %d\n" % (nspace))
    fout.write("veclen = %d\n" %(veclen))
    fout.write("data = %s\n" % (data_type))
    fout.write("field = %s\n" % (field_type))
    fout.write("label = %s\n" % (label))
    fout.write("\f\f")

    # write data
    print("num of elements = %d" % (len(pa)))
    for v in pa:
        value = float(v)
        fout.write(struct.pack("<f", value))
    
    # write coord
    print("num of coord = %d" % (len(coord)))
    for p in coord:
        value = float(p[0]) * 0.529177249
        fout.write(struct.pack("<f", value))
    for p in coord:
        value = float(p[1]) * 0.529177249
        fout.write(struct.pack("<f", value))
    for p in coord:
        value = float(p[2]) * 0.529177249
        fout.write(struct.pack("<f", value))

    fout.close()

if __name__ == '__main__':
    main()
