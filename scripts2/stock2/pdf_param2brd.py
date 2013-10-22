#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
transform ProteinDF parameter file to Bridge file.

This read the ProteinDF parameter file, which is usually named 'pdfparam.mpac',
and output Bridge-class coordination file formatted by MsgPack.
"""

import os
import optparse
import msgpack

from bratomgroup import *
from bratom import *
from brposition import *

def main():
    # parse options
    #
    parser = optparse.OptionParser(usage = "%prog [options]", version = "%prog 1.0")
    # set option
    parser.add_option("-w", "--write", dest="output_path",
                      help="output Bridge MsgPack file (default: bridge.mpac)", metavar="FILE",
                      default="bridge.mpac")
    parser.add_option("-v", "--verbose", dest="verbose",
                      help="verbose",
                      action="store_true", default=False)
    (opts, args) = parser.parse_args()
    
    pdfparam_path = "pdfparam.mpac"
    if (len(args) > 0):
        pdfparam_path = args.pop(0)
    brd_output_path = opts.output_path
    verbose = opts.verbose
    
    # read ProteinDF parameters
    if (verbose == True):
        print "loading: %s" % (pdfparam_path)
    f = open(pdfparam_path, "rb")
    contents = f.read()
    pdfparam_mpac = msgpack.unpackb(contents)
    f.close()

    # to Bridge class
    index = 0
    atom_group = BrAtomGroup()
    if (pdfparam_mpac.has_key("model") == True):
        model = pdfparam_mpac["model"]
        if (model.has_key("coordinates") == True):
            coord = model["coordinates"]
            for pdf_atom_group_key, pdf_atom_group in coord.items():
                subgroup = BrAtomGroup()
                for pdf_atom in pdf_atom_group:
                    element = pdf_atom.get("symbol", "")
                    position = pdf_atom.get("coord", [0.0, 0.0, 0.0])
                    charge = pdf_atom.get("charge", 0.0)
                    label = pdf_atom.get("label", None)
                    if (label == None):
                        label = ""
                    atom = BrAtom()
                    atom.set_element(element)
                    atom.set_position(BrPosition(position))
                    atom.set_name(label)
                    atom.set_charge(charge)
                    subgroup.set_atom(index, atom)
                    index += 1
                atom_group.set_group(pdf_atom_group_key, subgroup)

    # save
    bridge_mpac = msgpack.packb(atom_group.get_raw_data())
    savefile = open(brd_output_path, "wb")
    savefile.write(bridge_mpac)
    savefile.close()
            
        
if __name__ == '__main__':
    main()

        
