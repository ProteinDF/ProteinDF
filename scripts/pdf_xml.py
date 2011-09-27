#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import optparse
import re

from pdfdata import PdfData
from PdfOutput import PdfOutput
from PdfXml import PdfXml

pdf_work = "./fl_Work"

def main():
    # parse options
    #
    parser = optparse.OptionParser(usage="%prog ProteinDF_output", version="%prog 1.0")
    # set option
    parser.add_option("-v", action="store_true", dest="verbose", default=False)

    (opts, args) = parser.parse_args()

    if (len(args) <= 0):
        parser.print_help()
        exit()
    elif (len(args) >= 1):
        pdf_output_path = str(args[0])

    if (len(args) >= 2):
        pdf_work_path = str(args[1])

    if (opts.verbose == True):
        print "loading : %s" % (pdf_output_path)

    # set standard output data
    pdfOutput = PdfOutput(pdf_output_path)
    #print pdfOutput

    pdfdata = pdfOutput.get_pdfdata()
    #print pdfdata

    #pdf = PdfData()
    #pdf.setByPdfOutput(pdfOutput)

    # set energy level data
    method = pdfdata.get_method()
    set_energy_levels(method, pdfdata)

    # output XML
    pdfxml = PdfXml(pdfdata)
    xml = pdfxml.get_xml()
    print xml


def set_energy_levels(method, pdfdata):
    spin_list = set(['RKS'])
    if (method == 'sp'):
        spin_list = set(['UKS_ALPHA', 'UKS_BETA'])

    pdfconv = "%s/bin/pdfVector2txt" % (os.environ['PDF_HOME'])

    for iteration in range(1, (pdfdata.get_number_of_iterations() +1)):
        for spin in spin_list:
            suffix = ''
            if (spin == 'UKS_ALPHA'):
                suffix = 'a'
            elif (spin == 'UKS_BETA'):
                suffix = 'b'
            file_path = "%s/fl_Vct_Eigval%s%d" % (pdf_work, suffix, iteration)
            cmd = "%s -l %s 2>/dev/null" % (pdfconv, file_path)
            child = os.popen(cmd)
            data = child.read()
            err = child.close()
            #if (err):
            #    raise RuntimeError, '%r failed with exit code %d' % (cmd, err)

            data = re.sub("\n", ", ", str(data))
            data = re.sub(",\s*\n?$", "", data)
            pdfdata.set_energy_levels(iteration, spin, data)


if __name__ == '__main__':
    main()
    



