#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import optparse
import msgpack
import re

from dfdata import *

class DfDiffResults(object):
    def __init__(self, dfdata1, dfdata2, epsilon = 1.0E-5, verbose = False):
        assert(isinstance(dfdata1, DfData))
        assert(isinstance(dfdata2, DfData))
        self.dfdata1 = dfdata1
        self.dfdata2 = dfdata2
        self.isChecked = False
        self.report = ""
        self.isEquiv = True
        self.epsilon = epsilon
        self.verbose = verbose


    def check(self):
        if (self.isChecked == False):
            self.__check_number_of_atoms()
            self.__check_number_of_orbitals()
            self.__check_method()
            self.__check_occupation_level()
            self.__check_HOMO_level()
            self.__check_number_of_iterations()
            self.__check_total_energy()
            self.isChecked = True
        return self.isEquiv


    def __check_value(self, name, value1, value2):
        if (self.isEquiv == False):
            # skip check
            return
        
        if (isinstance(value1, str) == True):
            if (value1 != value2):
                self.isEquiv = False
                self.report += "[NG] (%s; \"%s\" != \"%s\")\n" % (name, value1, value2)
            elif (self.verbose == True):
                self.report += "[OK] (%s)\n" % (name)
            
        elif (isinstance(value1, float) == True):
            if (abs(value1 - value2) > self.epsilon):
                self.isEquiv = False
                self.report += "[NG] (%s; %d != %d)\n" % (name, value1, value2)
            elif (self.verbose == True):
                self.report += "[OK] (%s)\n" % (name)
        
        else:
            if (value1 != value2):
                self.isEquiv = False
                if (value1 == None):
                    value1 = "None"
                if (value2 == None):
                    value2 = "None"
                self.report += "[NG] (%s; %s != %s)\n" % (name, value1, value2)
            elif (self.verbose == True):
                self.report += "[OK] (%s)\n" % (name)
        

    def __check_number_of_atoms(self):
        self.__check_value("number of atoms",
                           self.dfdata1.get_number_of_atoms(),
                           self.dfdata2.get_number_of_atoms())

    def __check_number_of_orbitals(self):
        self.__check_value("number of orbitals",
                           self.dfdata1.get_number_of_orbitals(),
                           self.dfdata2.get_number_of_orbitals())

    def __check_method(self):
        self.__check_value("method",
                           self.dfdata1.get_method(),
                           self.dfdata2.get_method())

    def __check_occupation_level(self):
        method = self.dfdata1.get_method()
        if (method != self.dfdata2.get_method()):
            self.isEquiv = False
            return

        if (method == 'RKS'):
            self.__check_value("occupation_level(RKS)",
                               self.dfdata1.get_occupation_level('RKS'),
                               self.dfdata2.get_occupation_level('RKS'))
        elif ((method == 'UKS') or (method == 'ROKS')):
            self.__check_value("occupation_level(UKS_ALPHA)",
                               self.dfdata1.get_occupation_level('UKS_ALPHA'),
                               self.dfdata2.get_occupation_level('UKS_ALPHA'))
            self.__check_value("occupation_level(UKS_BETA)",
                               self.dfdata1.get_occupation_level('UKS_BETA'),
                               self.dfdata2.get_occupation_level('UKS_BETA'))
        else:
            self.isEquiv = False
            self.report = "NG [occupation_level] method type is wrong.\n"

        
    def __check_HOMO_level(self):
        method = self.dfdata1.get_method()
        if (method != self.dfdata2.get_method()):
            self.isEquiv = False
            return

        if (method == 'RKS'):
            self.__check_value("HOMO_level(RKS)",
                               self.dfdata1.get_HOMO_level('RKS'),
                               self.dfdata2.get_HOMO_level('RKS'))
        elif ((method == 'UKS') or (method == 'ROKS')):
            self.__check_value("HOMO_level(UKS_ALPHA)",
                               self.dfdata1.get_HOMO_level('UKS_ALPHA'),
                               self.dfdata2.get_HOMO_level('UKS_ALPHA'))
            self.__check_value("HOMO_level(UKS_BETA)",
                               self.dfdata1.get_HOMO_level('UKS_BETA'),
                               self.dfdata2.get_HOMO_level('UKS_BETA'))
        else:
            self.isEquiv = False
            self.report = "NG [HOMO_level] method type is wrong.\n"

        
    def __check_number_of_iterations(self):
        self.__check_value("number_of_iterations",
                           self.dfdata1.get_number_of_iterations(),
                           self.dfdata2.get_number_of_iterations())


    def __check_total_energy(self):
        max_itr = self.dfdata1.get_number_of_iterations()
        if (max_itr != self.dfdata2.get_number_of_iterations()):
            self.isEquiv = False
            return

        for itr in range(max_itr):
            name = "%d th total energy" % (itr +1)
            self.__check_value(name,
                               self.dfdata1.get_total_energy(itr),
                               self.dfdata2.get_total_energy(itr))
    

    def get_report(self):
        return self.report
        

def main():
    # parse options
    #
    parser = optparse.OptionParser(usage="%prog [options] mpac_file1 mpac_file2", version="%prog 1.0")
    # set option
    parser.add_option("-e", "--epsilon", dest="epsilon",
                      help="epsilon",
                      action="store", default=1.0E-5)
    parser.add_option("-v", "--verbose", dest="verbose",
                      help="verbose",
                      action="store_true")
    (opts, args) = parser.parse_args()

    # parameters
    if (len(args) < 2):
        parser.print_help()
        sys.exit(127)
        
    mpac_path1 = args[0]
    mpac_path2 = args[1]
    verbose = opts.verbose
    epsilon = float(opts.epsilon)

    # load & set
    dfdata1 = get_dfdata(mpac_path1, verbose)
    dfdata2 = get_dfdata(mpac_path2, verbose)

    # diff
    diff_obj = DfDiffResults(dfdata1, dfdata2, epsilon, verbose)
    isEquiv = diff_obj.check()

    sys.stderr.write(diff_obj.get_report())

    exit_value = 0
    if (isEquiv != True):
        sys.stderr.write("NOT equivalent\n")
        exit_value = 1
        
    sys.exit(exit_value)


def get_dfdata(mpac_path, verbose = False):
    if (verbose == True):
        sys.stderr.write("loading %s.\n" % (mpac_path))
    f = open(mpac_path, "rb")
    contents = f.read()
    raw_data = msgpack.unpackb(contents)
    f.close()
    
    dfdata = DfData()
    dfdata.set_raw_data(raw_data)

    return dfdata


if __name__ == '__main__':
    main()




