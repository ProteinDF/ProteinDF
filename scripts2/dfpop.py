#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import tempfile
import subprocess

from dfdata import *
from dfmatrix import *

class DfPopulation(object):
    """
    """
    __cmd = ""
    
    def __init__(self):
        """
        """
        self.__cmd = "%s/bin/mullikenPop" % (os.getenv('PDF_HOME'))

        
    def set_cmd(self, cmd):
        self.__cmd = cmd
        

    def get_atom_populations(self, iteration = 0):
        """
        """
        iteration = int(iteration)
        handle, filename = tempfile.mkstemp()
        os.close(handle)
        subprocess.check_call([self.__cmd, "-i", str(iteration), "-s", filename])
        
        mullikenPopMtx = DfMatrix()
        mullikenPopMtx.load_binary(filename)

        os.remove(filename)

        return mullikenPopMtx

        
if __name__ == '__main__':
    dfpop = DfPopulation()
    mtx = dfpop.get_atom_populations(0)
    mtx.printout()
    
