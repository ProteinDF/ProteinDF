#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import optparse

from PdfXml import PdfXml
from pdfdata import PdfData

def main():
    # parse options
    #
    parser = optparse.OptionParser(usage="%prog expected actual", version="%prog 1.0")
    # set option
    parser.add_option("-v", action="store_true", dest="verbose", default=False)

    (opts, args) = parser.parse_args()

    if (len(args) != 1):
        parser.print_help()
        sys.exit(1)

    xml_path = str(args[0])

    pdfxml = PdfXml(xml_path)
    #pdf = PdfData()
    #pdf.setByXmlFile(xmlPath)

    pdfdata = pdfxml.get_pdfdata()
    print pdfdata

    sys.exit(0)

def check(expect, actual, checkRange = 1.0E-5, verbose = False):
    errorCode = 0

    # model
    errorCode += checkInt(expect.getNumberOfAtoms(),
                          actual.getNumberOfAtoms(),
                          "number of atoms", verbose)
    errorCode += checkInt(expect.getNumberOfOrbitals(),
                          actual.getNumberOfOrbitals(),
                          "number of orbitals", verbose)
    
    # SCF
    errorCode += checkInt(expect.getNumberOfIterations(),
                          actual.getNumberOfIterations(),
                          "number of Iterations", verbose)
    errorCode += checkFloat(expect.getTotalEnergy(),
                            actual.getTotalEnergy(),
                            "total energy",
                            checkRange, verbose)
    if (expect.getNumberOfIterations() == actual.getNumberOfIterations()):
        numberOfIterations = expect.getNumberOfIterations()
        for i in range(1, numberOfIterations +1):
            errorCode += checkFloat(expect.getTotalEnergy(i),
                                    actual.getTotalEnergy(i),
                                    "%d th total energy" % (i),
                                    checkRange, verbose)

    return errorCode

def checkInt(expect, actual, descript, verbose = False):
    errorCode = 0
    judgeStr = ['pass', 'FAILED']

    if (expect != actual):
        errorCode = 1

    if ((errorCode != 0) or (verbose == True)):
        print "%s ... %s" % (descript, judgeStr[errorCode])

    return errorCode

def checkFloat(expect, actual, descript, checkRange = 1.0E-5, verbose = False):
    errorCode = 0
    judgeStr = ['pass', 'FAILED']

    if (abs(expect - actual) > checkRange):
        errorCode = 1

    if ((errorCode != 0) or (verbose == True)):
        print "%s ... %s" % (descript, judgeStr[errorCode])

    return errorCode


if __name__ == '__main__':
    main()
