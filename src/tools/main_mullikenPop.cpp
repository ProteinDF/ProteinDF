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

    std::string paramPath = "pdfparam.mpac";
    if (opt["p"].empty() != true) {
        paramPath = opt["p"];
    }
    const std::string savePath = opt["s"];

    const bool isVerbose = (opt["v"] == "defined");
    const bool showHelp = (opt["h"] == "defined");
    
    if (showHelp == true) {
        usage();
        return EXIT_SUCCESS;
    }


    TlMsgPack mpac;
    mpac.load(paramPath);
    TlSerializeData param = mpac.getSerializeData();

    int iteration = 0;
    if (opt["i"].empty() != true) {
        iteration = std::atoi(opt["i"].c_str());
    } else {
        iteration = param["iterations"].getInt();
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
    std::cerr <<  "  -s path:   save atom population data as pdf matrix file\n"; 
    std::cerr <<  "  -h:        show help(this)\n"; 
    std::cerr <<  "  -v:        verbose\n"; 
}

