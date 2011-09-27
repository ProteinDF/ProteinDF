#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
shrink ProteinDF output file.

This read the ProteinDF output file, which is usually named 'fl_Out_Std',
and output shrink data formatted by MsgPack.
"""

import optparse
import msgpack

from dfdata import *
from dfoutput import *
from dfvector import *

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
                      help="output MsgPack file (default: dfdata.mpac)", metavar="FILE",
                      default="dfdata.mpac")
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
    output = DfOutput(pdf_output_path)
    #print output

    # get DfData object
    dfdata = output.get_dfdata()
    #print dfdata

    # set energy level data
    method = dfdata.get_method()
    set_energy_levels(pdf_work_dir, method, dfdata)

    # output DfData as MsgPack
    mpac = msgpack.packb(dfdata.get_raw_data())

    # output file
    savefile = open(save_path, "wb");
    savefile.write(mpac)
    savefile.close()


def set_energy_levels(pdf_work_dir, method, dfdata):
    spin_list = set(['RKS'])
    if (method == 'sp'):
        spin_list = set(['UKS_ALPHA', 'UKS_BETA'])

    #pdfconv = "%s/bin/pdfVector2txt" % (os.environ['PDF_HOME'])

    for iteration in range(1, (dfdata.get_number_of_iterations() +1)):
        for spin in spin_list:
            suffix = ''
            if (spin == 'UKS_ALPHA'):
                suffix = 'a'
            elif (spin == 'UKS_BETA'):
                suffix = 'b'
            file_path = "%s/fl_Vct_Eigval%s%d" % (pdf_work_dir, suffix, iteration)
            #cmd = "%s -l %s 2>/dev/null" % (pdfconv, file_path)
            #child = os.popen(cmd)
            #data = child.read()
            #err = child.close()
            #if (err):
            #    raise RuntimeError, '%r failed with exit code %d' % (cmd, err)

            #data = re.sub("\n", ", ", str(data))
            #data = re.sub(",\s*\n?$", "", data)

            v = DfVector()
            v.load_binary(file_path)
            #dfdata.set_energy_levels(iteration, spin, v.vector.tolist())
            dfdata.set_energy_levels(iteration, spin, v.vector)

if __name__ == '__main__':
    main()
    



