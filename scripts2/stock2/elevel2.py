#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import getopt

from PdfXml import PdfXml
from pdfdata import PdfData
from graph_matplotlib import Graph

def main():
    # initialize
    output_file = "ks_level2"

    # parse args
    try:
        optlist, args = getopt.gnu_getopt(sys.argv[1:], "h",
                                          longopts=["help", "output="])
    except getopt.GetoptError:
        # exit before display the help message
        usage()
        sys.exit(1)
        
    for opt, arg in optlist:
        if opt in ("-h", "--help"):
            usage()
            sys.exit(0)
        if opt in ("--output"):
            output_file = arg
    
    if (len(args) == 0):
        usage()
        sys.exit(1)

    # setting
    xml_path = args.pop(0)
    pdfxml = PdfXml(xml_path)
    pdfdata = pdfxml.get_pdfdata()

    # plot
    graph = Graph(size=(8, 4))
    #graph = Graph()

    graph.set_xrange(-20.0, 5.0)
    graph.set_yrange(- 0.0, 2.0)
    graph.set_xlabel('energy level (eV)')
    #graph.draw_grid(False)

    HOMO_level = pdfdata.get_HOMO_level(method = 'RKS')
    energy_levels = pdfdata.get_energy_levels(
        spin_type = 'RKS',
        iteration = pdfdata.get_number_of_iterations())
    energy_level_list = energy_levels.split(",")
    # 1au = 27.21138 eV
    el = []
    for level in energy_level_list:
        el.append(float(level) * 27.21138)
        
    #print energy_level_list
    graph.plot_energy_levels_single(
        levels = el,
        homo_level = HOMO_level)

    graph.prepare()
    graph.file_out(output_file)


if __name__ == '__main__':
    main()




