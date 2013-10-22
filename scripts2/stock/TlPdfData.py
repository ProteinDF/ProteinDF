#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

class ProteinDF(object):
    """ProteinDF class"""
    # variables
    logLevel = 0
    numberOfAtoms = None
    numberOfOrbitals = None
    numberOfIterations = None
    scfData = {}

    # public
    def __init__(self, data = "", options = {}):
        self.logLevel = options.setdefault('log_level', 0)

    def getNumberOfAtoms(self):
        return self.numberOfAtoms

    def getNumberOfOrbitals(self):
        return self.numberOfOrbitals

    def setByPdfXml(self, xmlObject):
        self.numberOfAtoms = xmlObject.getNumberOfAtoms()
        self.numberOfOrbitals = xmlObject.getNumberOfOrbitals()

