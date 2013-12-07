#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import optparse
import msgpack
import re

from dfdata import *
from tlgraph import *

class PdfReport(object):
    def __init__(self, dfdata):
        self.dfdata = dfdata

    def getHtml(self):
        html = '''
<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<meta content="text/html; charset=utf-8" http-equiv="Content-Type" />
<title>ProteinDF report</title>
</head>

<body>


<table>
<tr>
<td width="400px">
<table style="width: 100%">

<tr>
<th>number of atoms</th>
<td>NUM_OF_ATOMS</td>
</tr>

<tr>
<th>number of orbitals</th>
<td>NUM_OF_AOS</td>
</tr>

<tr>
<th>Total Energy</th>
<td>TOTAL_ENERGY</td>
</tr>

<tr>
<th>HOMO-LUMO gap / eV</th>
<td>HOMO_LUMO_GAP</td>
</tr>

<tr>
<th>LUMO+4 / eV</th>
<td>LUMO_4_LEVEL</td>
</tr>

<tr>
<th>LUMO+3 / eV</th>
<td>LUMO_3_LEVEL</td>
</tr>

<tr>
<th>LUMO+2 / eV</th>
<td>LUMO_2_LEVEL</td>
</tr>

<tr>
<th>LUMO+1 / eV</th>
<td>LUMO_1_LEVEL</td>
</tr>

<tr>
<th>LUMO / eV</th>
<td>LUMO_0_LEVEL</td>
</tr>

<tr>
<th>HOMO / eV</th>
<td>HOMO_0_LEVEL</td>
</tr>

<tr>
<th>HOMO-1 / eV</th>
<td>HOMO_1_LEVEL</td>
</tr>

<tr>
<th>HOMO-2 / eV</th>
<td>HOMO_2_LEVEL</td>
</tr>

<tr>
<th>HOMO-3 / eV</th>
<td>HOMO_3_LEVEL</td>
</tr>

<tr>
<th>HOMO-4 / eV</th>
<td>HOMO_4_LEVEL</td>
</tr>

</table>
</td>
<td width="400px"><a href="./graph_elevel.png"><img src="./graph_elevel.png" style="width: 400px"/></a></td>
</tr>
<tr>
<td width="400px"><a href="./graph_TE.png"><img src="./graph_TE.png" style="width: 400px"/></a></td>
<td width="400px"><a href="./graph_convergence.png"><img src="./graph_convergence.png" style="width: 400px"/></a></td>
</tr>
</table>
</body>
</html>
'''
        # replace ==============================================================
        NUM_OF_ATOMS_str = "%d" % (self.dfdata.get_number_of_atoms())
        html = html.replace('NUM_OF_ATOMS', NUM_OF_ATOMS_str)

        NUM_OF_AOS_str = "%d" % (self.dfdata.get_number_of_orbitals())
        html = html.replace('NUM_OF_AOS', NUM_OF_AOS_str)

        total_energy = self.dfdata.get_total_energy()
        if (total_energy != None):
            TOTAL_ENERGY_str = "% 18.8f" % (total_energy)
            html = html.replace('TOTAL_ENERGY', TOTAL_ENERGY_str)
        else:
            html = html.replace('TOTAL_ENERGY', "N/A")

        # HOMO-LUMO gap
        HOMO_LUMO_gap = self.dfdata.get_band_gap()
        HOMO_LUMO_GAP_str = ""
        if (HOMO_LUMO_gap != None):
            HOMO_LUMO_GAP_str = "%2.3f" % (HOMO_LUMO_gap * 27.21138)
        html = html.replace('HOMO_LUMO_GAP', HOMO_LUMO_GAP_str)

        # 
        HOMO_level = self.dfdata.get_HOMO_level(method = 'RKS')
        levels = self.dfdata.get_energy_levels()
        for i in range(5):
            replace_str = "HOMO_%d_LEVEL" % (i)
            if (((HOMO_level - i) <= len(levels)) and
                ((HOMO_level - i) >= 0)):
                html = html.replace(replace_str, "%2.3f" % (levels[HOMO_level -1 -i] * 27.21138)) # option base 0
            else:
                html = html.replace(replace_str, "N/A")

        for i in range(5):
            replace_str = "LUMO_%d_LEVEL" % (i)
            if (((HOMO_level + 1 +i) <= len(levels)) and
                ((HOMO_level + 1 +i) >= 0)):
                html = html.replace(replace_str, "%2.3f" % (levels[HOMO_level -1 +1 +i] * 27.21138)) # option base 0
            else:
                html = html.replace(replace_str, "N/A")

        return html

def main():
    # parse options
    #
    parser = optparse.OptionParser(usage="%prog [options] <mpac file>", version="%prog 1.0")
    # set option
    parser.add_option("-f", "--file", dest="graph_path",
                      help="output graph base file name (default: graph)", metavar="FILE",
                      default="graph")
    parser.add_option("-v", "--verbose", dest="verbose",
                      help="verbose",
                      action="store_true", default=False)
    (opts, args) = parser.parse_args()

    # variable
    mpac_path = "results.mpac"
    if (len(args) > 0):
        mpac_path = args.pop(0)
    graph_path = opts.graph_path
    ext = "png" # graph file extension
    verbose = opts.verbose

    # load
    if (verbose == True):
        sys.stderr.write("loading: %s.\n" % (mpac_path))
    f = open(mpac_path, "rb")
    contents = f.read()
    raw_data = msgpack.unpackb(contents)
    f.close()
    
    dfdata = DfData()
    dfdata.set_raw_data(raw_data)

    # plot (convergence MO levels)
    conv_elevel_path = "%s_elevel.%s" % (graph_path, ext)
    if (verbose == True):
        sys.stderr.write("plot convergence MO levels: %s.\n" % (conv_elevel_path));
    plot_convergence_energy_level(dfdata, conv_elevel_path)

    # plot (convergence last MO levels)
    last_elevel_path = "%s_elevel_last.%s" % (graph_path, ext)
    if (verbose == True):
        sys.stderr.write("plot last energy level: %s.\n" % (last_elevel_path));
    plot_energy_level(dfdata, last_elevel_path)

    # plot (convergence history)
    convergence_history_path = "%s_convergence.%s" % (graph_path, ext)
    if (verbose == True):
        sys.stderr.write("plot convergence history: %s.\n" % (convergence_history_path));
    plot_convergence_check(dfdata, convergence_history_path)

    # plot (total energy)
    convergence_TE_path = "%s_TE.%s" % (graph_path, ext)
    if (verbose == True):
        sys.stderr.write("plot total energy: %s.\n" % (convergence_TE_path));
    plot_convergence_TE(dfdata, convergence_TE_path)


    # HTML
    report = PdfReport(dfdata)
    html_file = open('report.html', 'w')
    html_file.write(report.getHtml())
    html_file.close()

    sys.exit()
    

def plot_convergence_energy_level(dfdata, output_path):
    graph = TlGraphEnergyLevelHistory()

    graph.set_xrange(0, dfdata.get_number_of_iterations())
    graph.set_yrange(-20.0, 5.0)
    graph.set_xlabel('SCF convergence step')
    graph.set_ylabel('energy level / eV')
    graph.draw_grid(True)

    HOMO_level = dfdata.get_HOMO_level(method = 'RKS')
    for itr in range(1, dfdata.get_number_of_iterations() +1):
        energy_levels = dfdata.get_energy_levels(
            spin_type = 'RKS',
            iteration = itr,
            request_unit = "eV")
        
        graph.plot_energy_levels(iteration = itr,
                                 levels = energy_levels,
                                 homo_level = HOMO_level -1)

    graph.prepare()
    graph.file_out(output_path)


def plot_energy_level(dfdata, output_path):
    graph = TlGraphEnergyLevelSingle(size=(8, 4))

    graph.set_xrange(-20.0, 5.0)
    graph.set_yrange(- 0.0, 2.0)
    graph.set_xlabel('energy level / eV')
    
    HOMO_level = dfdata.get_HOMO_level(method = 'RKS')
    energy_levels = dfdata.get_energy_levels(
        spin_type = 'RKS',
        iteration = dfdata.get_number_of_iterations(),
        request_unit = "eV")

    graph.plot_energy_levels(levels = energy_levels,
                             homo_level = HOMO_level -1)
    
    graph.prepare()
    graph.file_out(output_path)


def plot_convergence_check(dfdata, output_path):
    graph = TlGraphConvergenceCheck()

    graph.set_xlabel('SCF convergence step')
    graph.set_ylabel('difference / a.u.')
    graph.draw_grid(True)
    graph.draw_legend(True)

    iterations = dfdata.get_number_of_iterations()
    # total energy
    TE_history = range(iterations +1)
    TE_history[0] = None
    # delta TE
    deltaTE_history = range(iterations +1)
    deltaTE_history[0] = None
    # delta Density Matrix
    deltaDensMat_history = range(iterations +1)
    deltaDensMat_history[0] = None
    for itr in range(1, iterations +1):
        TE_history[itr] = dfdata.get_total_energy(itr)
        info = dfdata.get_convergence_info(itr)
        deltaTE_history[itr] = info.get('max_deviation_of_total_energy', None)
        deltaDensMat_history[itr] = info.get('max_deviation_of_density_matrix', None)

    #graph.plot(TE_history, "Total Energy")
    graph.plotLog(deltaTE_history, "delta Total Energy")
    graph.plotLog(deltaDensMat_history, "delta Density Matrix")

    graph.prepare()
    graph.file_out(output_path)


def plot_convergence_TE(dfdata, output_path):
    graph = TlGraphConvergenceCheck()

    graph.set_xlabel('SCF convergence step')
    graph.set_ylabel('Total Energy / a.u.')
    graph.draw_grid(True)
    graph.draw_legend(True)

    iterations = dfdata.get_number_of_iterations()
    # total energy
    TE_history = range(iterations +1)
    TE_history[0] = None
    for itr in range(1, iterations +1):
        TE_history[itr] = dfdata.get_total_energy(itr)

    graph.plot(TE_history, "Total Energy")

    graph.prepare()
    graph.file_out(output_path)


if __name__ == '__main__':
    main()




