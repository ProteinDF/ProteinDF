#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
extract ProteinDF output files.

This read the ProteinDF output file, which is usually named 'fl_Out_Std',
and output shrink data formatted by MsgPack.
"""

import os
import optparse
import msgpack

from dfdata import *
from dfoutput import *
from dfvector import *
from dfmatrix import *

def main():
    # parse options
    #
    parser = optparse.OptionParser(usage = "%prog [options]", version = "%prog 1.0")
    # set option
    parser.add_option("-f", "--file", dest="dfoutput_path",
                      help="ProteinDF output file (default: fl_Out_Std)", metavar="FILE",
                      default="fl_Out_Std")
    parser.add_option("-d", "--directory", dest="directory",
                      help="ProteinDF calculation directory (default: current)", metavar="DIR",
                      default=".")
    parser.add_option("-w", "--write", dest="mpac_path",
                      help="output MsgPack file (default: results.mpac)", metavar="FILE",
                      default="results.mpac")
    parser.add_option("-v", "--verbose", dest="verbose",
                      help="verbose",
                      action="store_true", default=False)
    (opts, args) = parser.parse_args()

    pdf_output_path = opts.dfoutput_path
    pdf_calc_dir = opts.directory
    pdf_work_dir = "%s/fl_Work" % (pdf_calc_dir)
    save_path = opts.mpac_path

    if (opts.verbose == True):
        print "loading : %s" % (pdf_output_path)

    # set standard output data
    df_output = DfOutput(pdf_output_path)
    dfdata = df_output.get_dfdata()
    #print dfdata

    # read ProteinDF input parameters
    pdfparam_path = "pdfparam.mpac"
    dfdata_input = get_dfdata_from_dfparam(pdfparam_path)
    dfdata.update(dfdata_input)

    # set energy level data
    method = dfdata.get_method()
    set_energy_levels(pdf_work_dir, method, dfdata)

    # set LCAO of last iteration
    set_lcao_matrix(pdf_work_dir, method, dfdata)

    # output DfData as MsgPack
    mpac = msgpack.packb(dfdata.get_raw_data())

    # output file
    savefile = open(save_path, "wb");
    savefile.write(mpac)
    savefile.close()


# ==============================================================================
def get_dfdata_from_dfparam(dfparam_path):
    dfdata = DfData()
    # loading
    if (os.path.isfile(dfparam_path) == True):
        f = open(dfparam_path, "rb")
        contents = f.read()
        data = msgpack.unpackb(contents)
        f.close()
        
        if (data.has_key("model") == True):
            model = data["model"]
            if (model.has_key("atoms") == True):
                dfdata.set_number_of_atoms(model["atoms"])
            if (model.has_key("AOs") == True):
                dfdata.set_number_of_orbitals(model["AOs"])
            if (model.has_key("method") == True):
                method = model["method"]
                dfdata.set_method(method)
                if (dfdata.get_method() == "RKS"):
                    if (model.has_key("method/nsp/occlevel") == True):
                        dfdata.set_occupation_level("RKS", model["method/nsp/occlevel"])
                elif (dfdata.get_method() == "UKS"):
                    if (model.has_key("method/sp/occlevel") == True):
                        dfdata.set_occupation_level("UKS_ALPHA", model["method/sp/alpha-spin-occlevel"])
                    if (model.has_key("method/sp/occlevel") == True):
                        dfdata.set_occupation_level("UKS_BETA", model["method/sp/beta-spin-occlevel"])
                else:
                    print("unsupport method: %s" % (method))
            if (model.has_key("iterations") == True):
                dfdata.set_number_of_iterations(model["iterations"]);

    return dfdata


def set_energy_levels(pdf_work_dir, method, dfdata):
    spin_list = set()
    if (method == 'RKS'):
        spin_list = set(['RKS'])
    elif ((method == 'UKS') or (method == 'ROKS')):
        spin_list = set(['UKS_ALPHA', 'UKS_BETA'])

    for iteration in range(1, (dfdata.get_number_of_iterations() +1)):
        for spin in spin_list:
            suffix = ''
            if (spin == 'UKS_ALPHA'):
                suffix = 'a'
            elif (spin == 'UKS_BETA'):
                suffix = 'b'
            file_path = "%s/fl_Vct_Eigval%s%d" % (pdf_work_dir, suffix, iteration)

            v = DfVector()
            v.load_binary(file_path)
            dfdata.set_energy_levels(iteration, spin, v.vector)


def set_lcao_matrix(pdf_work_dir, method, dfdata):
    spin_list = set()
    if (method == 'RKS'):
        spin_list = set(['RKS'])
    elif ((method == 'UKS') or (method == 'ROKS')):
        spin_list = set(['UKS_ALPHA', 'UKS_BETA'])

    last_itr = dfdata.get_number_of_iterations()
    for spin in spin_list:
        if (spin == 'RKS'):
            suffix = 'rks'
        if (spin == 'UKS_ALPHA'):
            suffix = 'uks-alpha'
        elif (spin == 'UKS_BETA'):
            suffix = 'uks-beta'
    file_path = "%s/fl_Mtr_C.matrix.%s%d" % (pdf_work_dir, suffix, last_itr)

    m = DfMatrix()
    m.load_binary(file_path)
    dfdata.set_lcao_matrix(last_itr, spin, m)


if __name__ == '__main__':
    main()
    
