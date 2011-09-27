#include <cassert>
#include <iostream>
#include <cstdlib>
#include <string>

#include "DfPopulation.h"
#include "TlGetopt.h"
#include "TlVector.h"
#include "TlSymmetricMatrix.h"
#include "TlMsgPack.h"
#include "TlSerializeData.h"

void usage();

int main(int argc, char* argv[])
{
    TlGetopt opt(argc, argv, "hi:vs:");
    const bool showHelp = (opt["h"] == "defined");
    const bool isVerbose = (opt["v"] == "defined");
    const std::string savePath = opt["s"];
    
    if (showHelp == true) {
        usage();
        return EXIT_SUCCESS;
    }
    
    TlMsgPack mpac;
    mpac.load("pdfparam.mpac");
    TlSerializeData param = mpac.getSerializeData();

    int iteration = 0;
    if (opt["i"].empty() != true) {
        iteration = std::atoi(opt["i"].c_str());
    } else {
        iteration = param["model"]["iterations"].getInt();
    }

    if (isVerbose == true) {
        std::cerr << "iteration = " << iteration << std::endl;
    }
    
    DfPopulation dfPop(&param);

    if (savePath.empty() != true) {
        const TlMatrix mtx = dfPop.getAtomPopData(iteration);
        if (isVerbose == true) {
            std::cerr << "save Mulliken Data as " << savePath << std::endl;
        }
        mtx.save(savePath);
    } else {
        dfPop.getReport(iteration, std::cout);
    }

    return EXIT_SUCCESS;
}


void usage()
{
    std::cerr << "calculation Mulliken Population." << std::endl;
    std::cerr << "usage: mullikenPop [options] " << std::endl;
    std::cerr <<  "  -i num:    set SCF iteration to get Mulliken Population\n"; 
    std::cerr <<  "  -s path:   save population data file\n"; 
    std::cerr <<  "  -h:        show help(this)\n"; 
    std::cerr <<  "  -v:        verbose\n"; 
}

