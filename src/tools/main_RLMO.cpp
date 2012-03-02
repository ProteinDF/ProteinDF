#include <iostream>
#include <cstdlib>

#include "DfRLMO.h"
#include "TlGetopt.h"
#include "TlSerializeData.h"
#include "TlMsgPack.h"
#include "TlOrbitalInfo.h"

void showHelp()
{
    std::cout << "RLMO [options] start_block_AO_index ..." << std::endl;
    std::cout << " OPTIONS:" << std::endl;
    std::cout << "  -h:      show help" << std::endl;
    std::cout << "  -v:      verbose" << std::endl;
}

int main(int argc, char* argv[])
{
    TlGetopt opt(argc, argv, "hv");
    
    if (opt["h"] == "defined") {
        showHelp();
        return EXIT_SUCCESS;
    }
    
    const bool bVerbose = (opt["v"] == "defined");

    TlMsgPack mpac;
    mpac.load("pdfparam.mpac");
    TlSerializeData param = mpac.getSerializeData();

    std::vector<DfObject::index_type> startBlockAOs;
    if (opt.getCount() > 1) {
        int numOfBlocks = opt.getCount() -1;
        for (int i = 0; i < numOfBlocks; ++i) {
            startBlockAOs.push_back(std::atoi(opt[i +1].c_str()));
        }
    } else {
        TlOrbitalInfo orbInfo(param["coordinates"],
                              param["basis_sets"]);
        int atomIndex = -1;
        const std::size_t numOfAOs = orbInfo.getNumOfOrbitals();
        for (DfObject::index_type i = 0; i < numOfAOs; ++i) {
            int currentAtomIndex = orbInfo.getAtomIndex(i);
            if (atomIndex != currentAtomIndex) {
                startBlockAOs.push_back(i);
                atomIndex = currentAtomIndex;
            }
        }
    }
    std::cerr << "number of blocks = " << startBlockAOs.size() << std::endl;
    
    DfRLMO dfRLMO(&param);
    dfRLMO.exec(startBlockAOs);

    return EXIT_SUCCESS;
}


