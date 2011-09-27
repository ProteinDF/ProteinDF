#!/usr/bin/env python
# -*- coding: utf-8 -*-

#from pylab import figure, show
from Numeric import *
from matplotlib.ticker import NullLocator

#--- for non-X environment ---
# set  matplotlib display before "import pylab"
import matplotlib
matplotlib.use('Agg')
import pylab

class Graph(object):
    __param = {}
    __fig = None
    __ax = None

    def __init__(self, size = None):
        self.__fig = pylab.figure(figsize = size)
        self.__ax = self.__fig.add_subplot(111)
        self.__title = ""

        # clear the figure
        #self.__fig.clf()
    
    def file_out(self, file_path):
        #self.prepare()
        self.__fig.savefig(file_path)


    def set_xrange(self, min_value, max_value):
        self.__param['min_x'] = min_value
        self.__param['max_x'] = max_value


    def set_yrange(self, min_value, max_value):
        self.__param['min_y'] = min_value
        self.__param['max_y'] = max_value


    def set_xlabel(self, label):
        self.__param['xlabel'] = label


    def set_ylabel(self, label):
        self.__param['ylabel'] = label


    def draw_grid(self, yn):
        assert((yn == True) or (yn == False))
        self.__param['draw_grid'] = yn


    def plot_energy_levels(self, iteration, levels, homo_level =-1):
        """
        levels is the list of energy level.
        """
        level = 1
        for level_value in levels:
            level_value = float(level_value)
            if (level != homo_level):
                self.__ax.broken_barh([(iteration -0.45, 0.9)], (level_value, 0),
                                      edgecolors='black', facecolors='black')
                                      #linewidth='1')
            else:
                self.__ax.broken_barh([(iteration -0.50, 1.0)], (level_value, 0),
                                      edgecolors='red', facecolors='red')
                                      #linewidth='1')
            level += 1

    def plot_energy_levels_single(self, levels, homo_level =-1):
        """
        levels is the list of energy level.
        """
        level = 1
        self.__param['ylocator'] = False
        for level_value in levels:
            level_value = float(level_value)
            if (level != homo_level):
                self.__ax.broken_barh([(level_value, 0)], (1.0 -0.50, 1.0),
                                      edgecolors='black', facecolors='black'
                                      )
            else:
                self.__ax.broken_barh([(level_value, 0)], (1.0 -0.50, 1.0),
                                      edgecolors='red', facecolors='red'
                                      )
            level += 1


    def plot_basis(self, type = 's', coef = [], exponent = [], color = 'blue'):
        """
        print Gaussian
        """
        assert(len(coef) == len(exponent))

        max_value = 5
        x = []
        for i in range(max_value * 100):
            v = 0.01 * float(i)
            x.append(v)

        y = []
        if (type == 's'):
            for i in x:
                i = float(i)
                numOfItems = len(coef)
                v = 0.0
                for t in range(numOfItems):
                    tmp = float(exponent[t]) * (i * i)
                    if (fabs(tmp) > log(10000.0)):
                        continue
                    tmp2 = exp(- tmp)
                    v += float(coef[t]) * tmp2
                y.append(v)
        elif (type == 'p'):
            for i in x:
                i = float(i)
                numOfItems = len(coef)
                v = 0.0
                for t in range(numOfItems):
                    tmp = float(exponent[t]) * (i * i)
                    if (fabs(tmp) > log(10000.0)):
                        continue
                    tmp2 = exp(- tmp)
                    v += float(coef[t]) * tmp2 * i
                y.append(v)
        elif (type == 'd'):
            for i in x:
                i = float(i)
                numOfItems = len(coef)
                v = 0.0
                for t in range(numOfItems):
                    tmp = float(exponent[t]) * (i * i)
                    if (fabs(tmp) > log(10000.0)):
                        continue
                    tmp2 = exp(- tmp)
                    v += float(coef[t]) * tmp2 * i * i
                y.append(v)

        self.__ax.plot(x, y, color)


    def prepare(self):
        # range
        if (self.__param.has_key('max_x') == True):
            self.__ax.set_xlim(self.__ax.get_xlim()[0], self.__param['max_x'])
        if (self.__param.has_key('min_x') == True):
            self.__ax.set_xlim(self.__param['min_x'], self.__ax.get_xlim()[1])
        if (self.__param.has_key('max_y') == True):
            self.__ax.set_ylim(self.__ax.get_ylim()[0], self.__param['max_y'])
        if (self.__param.has_key('min_y') == True):
            self.__ax.set_ylim(self.__param['min_y'], self.__ax.get_ylim()[1])

        if (self.__param.get('ylocator', True) == False):
            nullLocator = NullLocator()
            self.__ax.yaxis.set_major_locator(nullLocator)
            self.__ax.yaxis.set_minor_locator(nullLocator)

        # axis label
        if (self.__param.has_key('xlabel') == True):
            self.__ax.set_xlabel(self.__param['xlabel'])

        if (self.__param.has_key('ylabel') == True):
            self.__ax.set_ylabel(self.__param['ylabel'])

        # grid
        if (self.__param.get('draw_grid', False) == True):
            self.__ax.grid(True)

