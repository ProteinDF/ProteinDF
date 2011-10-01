#!/usr/bin/env python
# -*- coding: utf-8 -*-

#from pylab import figure, show
import pylab
from Numeric import *
from matplotlib.ticker import NullLocator

class TlGraph(object):
    def __init__(self, size = None):
        self.param = {}
        self.fig = pylab.figure(figsize = size)
        self.ax = self.fig.add_subplot(111)
        self.__title = ""

        # clear the figure
        #self.fig.clf()
    
    def file_out(self, file_path):
        #self.prepare()
        self.fig.savefig(file_path)


    def set_xrange(self, min_value, max_value):
        self.param['min_x'] = min_value
        self.param['max_x'] = max_value


    def set_yrange(self, min_value, max_value):
        self.param['min_y'] = min_value
        self.param['max_y'] = max_value


    def set_xscale(self, value):
        assert(value in ['linear', 'log', 'symlog'])
        self.param['xscale'] = value
 

    def set_yscale(self, value):
        assert(value in ['linear', 'log', 'symlog'])
        self.param['yscale'] = value
        

    def set_xlabel(self, label):
        self.param['xlabel'] = label


    def set_ylabel(self, label):
        self.param['ylabel'] = label


    def draw_grid(self, yn):
        assert((yn == True) or (yn == False))
        self.param['draw_grid'] = yn


    def draw_legend(self, yn):
        assert((yn == True) or (yn == False))
        self.param['draw_legend'] = yn


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

        self.ax.plot(x, y, color)


    def prepare(self):
        # range
        if (self.param.has_key('max_x') == True):
            self.ax.set_xlim(self.ax.get_xlim()[0], self.param['max_x'])
        if (self.param.has_key('min_x') == True):
            self.ax.set_xlim(self.param['min_x'], self.ax.get_xlim()[1])
        if (self.param.has_key('max_y') == True):
            self.ax.set_ylim(self.ax.get_ylim()[0], self.param['max_y'])
        if (self.param.has_key('min_y') == True):
            self.ax.set_ylim(self.param['min_y'], self.ax.get_ylim()[1])

        # scale
        if (self.param.has_key('xscale') == True):
            self.ax.set_xscale(self.param['xscale'])
        if (self.param.has_key('yscale') == True):
            self.ax.set_yscale(self.param['yscale'])

        # locator
        if (self.param.get('xlocator', True) == False):
            nullLocator = NullLocator()
            self.ax.xaxis.set_major_locator(nullLocator)
            self.ax.xaxis.set_minor_locator(nullLocator)
        if (self.param.get('ylocator', True) == False):
            nullLocator = NullLocator()
            self.ax.yaxis.set_major_locator(nullLocator)
            self.ax.yaxis.set_minor_locator(nullLocator)

        # axis label
        if (self.param.has_key('xlabel') == True):
            self.ax.set_xlabel(self.param['xlabel'])
        if (self.param.has_key('ylabel') == True):
            self.ax.set_ylabel(self.param['ylabel'])

        # grid
        if (self.param.setdefault('draw_grid', False) == True):
            self.ax.grid(True)

        # legend
        if (self.param.setdefault('draw_legend', False) == True):
            self.ax.legend()


class TlGraphEnergyLevelHistory(TlGraph):
    def __init__(self, size = None):
        TlGraph.__init__(self, size)

    def plot_energy_levels(self, iteration, levels, homo_level = -1):
        """
        levels is the list of energy level.
        """
        if (levels == None):
            return
        
        min_y = self.param['min_y']
        max_y = self.param['max_y']

        # this value is for negligence
        last_value = min_y
        negligence_range = (max_y - min_y) * 0.001 # 0.1%
        
        for level, value in enumerate(levels):
            # draw HOMO
            if (level == homo_level):
                self.ax.broken_barh([(iteration -0.50, 1.0)], (value, 0),
                                    edgecolors='red', facecolors='red')
                continue
            
            # negligence
            if ((value < min_y) or (max_y < value)):
                continue
            if (abs(last_value - value) < negligence_range):
                continue
            else:
                last_value = value

            # register data
            self.ax.broken_barh([(iteration -0.45, 0.9)], (value, 0),
                                edgecolors='black', facecolors='black')


class TlGraphEnergyLevelSingle(TlGraph):
    def __init__(self, size = None):
        TlGraph.__init__(self, size)

    def plot_energy_levels(self, levels, homo_level = -1):
        """
        levels is the list of energy level.
        """
        self.param['ylocator'] = False

        min_y = self.param['min_y']
        max_y = self.param['max_y']

        # this value is for negligence
        last_value = min_y
        negligence_range = (max_y - min_y) * 0.001 # 0.1%
        
        # plot
        for level, value in enumerate(levels):
            # draw HOMO
            if (level == homo_level):
                self.ax.broken_barh([(value, 0)], (1.0 -0.50, 1.0),
                                    edgecolors='red', facecolors='red'
                                    )

            # negligence
            if ((value < min_y) or (max_y < value)):
                continue
            if (abs(last_value - value) < negligence_range):
                continue
            else:
                last_value = value

            # register data
            self.ax.broken_barh([(value, 0)], (1.0 -0.50, 1.0),
                                edgecolors='black', facecolors='black'
                                )


class TlGraphTotalEnergyHistory(TlGraph):
    def __init__(self, size = None):
        TlGraph.__init__(self, size)

    def plot(self, values, label = "", marker = '+', linestyle = '-'):
        """
        plot total energy history data
        """
        x = []
        y = []
        for i in range(len(values)):
            if (values[i] != None):
                x.append(i)
                y.append(values[i])
        self.ax.plot(x, y, label = label, marker = marker, linestyle = linestyle)


class TlGraphConvergenceCheck(TlGraph):
    def __init__(self, size = None):
        TlGraph.__init__(self, size)

    def plot(self, values, label = "", marker = '+', linestyle = '-'):
        """
        plot convergence history data
        """
        x = []
        y = []
        for i in range(len(values)):
            if (values[i] != None):
                x.append(i)
                y.append(values[i])
        self.ax.plot(x, y, label = label, marker = marker, linestyle = linestyle)

    def plotLog(self, values, label = "", marker = '+', linestyle = '-'):
        """
        plot convergence history data (log)
        """
        x = []
        y = []
        for i in range(len(values)):
            if (values[i] != None):
                x.append(i)
                y.append(values[i])
        self.ax.semilogy(x, y, label = label, marker = marker, linestyle = linestyle)
