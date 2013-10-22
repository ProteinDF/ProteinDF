#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import optparse
import msgpack

from dfdata import *
from tlgraph import *

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
    mpac_path = "dfdata.mpac"
    if (len(args) > 0):
        mpac_path = args.pop(0)
    graph_path = opts.graph_path
    ext = "png" # graph file extension
    verbose = opts.verbose

    # load
    if (verbose == True):
        sys.stderr.write("loading %s.\n" % (mpac_path))
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
            unit = "a.u.")
        
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
        unit = "a.u.")

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




