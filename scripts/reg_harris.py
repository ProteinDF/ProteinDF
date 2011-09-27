#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
make msgpack files from the results of ProteinDF calculation for Harris functional DB.
"""

import os
import sys
import optparse
import msgpack
import array

import dfmatrix

def main():
    # initialize
    parser = optparse.OptionParser(usage="%prog FILE", version="%prog 1.0")
    parser.add_option("-o", "--output", dest="db_path", metavar="PATH",
                      default="harris.mpac",
                      help="update database file path")
    parser.add_option("-p", "--param", dest="param", metavar="PATH",
                      default="pdfparam.mpac",
                      help="ProteinDF parameter path")

    (opts, args) = parser.parse_args()

    # set up parameters
    matrixPathList = [] # density matrix list
    if (len(args) > 0):
        for arg in args:
            matrixPathList.append(arg)
    else:
        parser.print_help()
        sys.exit()

    pdfparam = {}
    if (os.path.exists(opts.param)):
        print("ProteinDF param path: %s" % (opts.param))
        f = open(opts.param, "rb")
        contents = f.read()
        pdfparam = msgpack.unpackb(contents) # ProteinDF parameters
        f.close()
    else:
        print("could not open: %s" % (opts.param))
        parser.print_help()
        sys.exit()

    db = {}
    if (os.path.exists(opts.db_path)):
        print("load database: %s" % (opts.db_path))
        db_file = open(opts.db_path, "rb")
        db_contents = db_file.read()
        db = msgpack.unpackb(db_contents)
        db_file.close()

    # set basis_set
    db.setdefault('model', {})
    db['model'].setdefault('basis_set', {})
    atom = str(pdfparam['model']['basis_set'].keys().pop())
    db['model']['basis_set'][atom] = pdfparam['model']['basis_set'][atom]

    # coord
    db.setdefault('harris', {})
    db['harris'].setdefault(atom, {})
    db['harris'][atom] = pdfparam['model']['coordinates']

    # density matrix
    numOfAOs = pdfparam['model']['AOs']
    mat = dfmatrix.DfSymmetricMatrix(numOfAOs)
    for matrixPath in matrixPathList:
        print("add density matrix: %s" % (matrixPath))
        tmp = dfmatrix.DfSymmetricMatrix()
        tmp.load_binary(matrixPath)
        mat = mat + tmp
    db.setdefault('density_matrix', {})
    db['density_matrix'].setdefault(atom, {})
    db['density_matrix'][atom]['row'] = int(mat.get_num_of_rows())
    db['density_matrix'][atom]['col'] = int(mat.get_num_of_cols())
    db['density_matrix'][atom]['type'] = "RLHD"
    matrix_raw_data = mat.get_raw_data()
    db['density_matrix'][atom]['data'] = matrix_raw_data["data"]
    
    # output
    fout = open(opts.db_path, "wb")
    db_mpac = msgpack.packb(db)
    fout.write(db_mpac)
    fout.close()


if __name__ == '__main__':
    main()
