#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import optparse

from pdfdata import PdfData
from PdfXml import PdfXml

def main():
    # parse options
    #
    parser = optparse.OptionParser(usage="%prog expected actual", version="%prog 1.0")
    # set option
    parser.add_option("-v", action="store_true", dest="verbose", default=False)

    (opts, args) = parser.parse_args()

    if (len(args) != 2):
        parser.print_help()
        exit()

    expectXmlPath = str(args[0])
    actualXmlPath = str(args[1])

    expect_xml = PdfXml(expectXmlPath)
    actual_xml = PdfXml(actualXmlPath)

    #print expect_xml
    #print actual_xml

    expect = expect_xml.get_pdfdata()
    actual = actual_xml.get_pdfdata()

    #print expect
    #print ">>>>"
    #print actual

    # check
    checkRange = 1.0E-5
    errorCode = check(expect, actual, checkRange, opts.verbose)

    if (errorCode == 0):
        if (opts.verbose == True):
            print "normal terminate."
    else:
        print "error found."
        
    sys.exit(errorCode)


def check(expect, actual, checkRange = 1.0E-5, verbose = False):
    assert(isinstance(expect, PdfData))
    assert(isinstance(actual, PdfData))

    errorCode = 0

    # model
    errorCode += checkInt(expect.get_number_of_atoms(),
                          actual.get_number_of_atoms(),
                          "number of atoms", verbose)
    errorCode += checkInt(expect.get_number_of_orbitals(),
                          actual.get_number_of_orbitals(),
                          "number of orbitals", verbose)
    
    # SCF
    errorCode += checkInt(expect.get_number_of_iterations(),
                          actual.get_number_of_iterations(),
                          "number of Iterations", verbose)
    if (expect.get_number_of_iterations() == actual.get_number_of_iterations()):
        numberOfIterations = expect.get_number_of_iterations()
        for i in range(1, numberOfIterations +1):
            errorCode += checkFloat(expect.get_total_energy(i),
                                    actual.get_total_energy(i),
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
