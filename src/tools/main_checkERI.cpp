#include <iostream>
#include <cstdlib>
#include <vector>

#include "DfEriEngine.h"
#include "TlGetopt.h"
#include "TlUtils.h"
#include "TlFile.h"
#include "TlOrbitalInfo.h"
#include "TlSymmetricMatrix.h"
#include "TlMsgPack.h"
#include "TlSerializeData.h"
#include "TlCombineDensityMatrix.h"

void showHelp()
{
    std::cout << "combineDensityMatrix PDF_REFERENCE_DIR PDF_TARGET_DIR [UPDATE_MATRIX_FILE]" << std::endl;
    std::cout << "OPTIONS:" << std::endl;
    std::cout << "  -c: check matrix" << std::endl;
    std::cout << "  -r <range>: range" << std::endl;
    std::cout << "  -x: disable position check" << std::endl;
    std::cout << "  -h: show help(this)" << std::endl;
    std::cout << "  -v: verbose" << std::endl;
}


int main(int argc, char* argv[])
{
    TlGetopt opt(argc, argv, "hv");

    // parameters
    const bool isVerbose = (opt["v"] == "defined");

    int shellIndexP = std::atoi(opt[1].c_str());
    int shellIndexQ = std::atoi(opt[2].c_str());
    int shellIndexR = std::atoi(opt[3].c_str());
    int shellIndexS = std::atoi(opt[4].c_str());
    std::cerr << TlUtils::format("(%d %d|%d %d)",
                                 shellIndexP, shellIndexQ, shellIndexR, shellIndexS)
              << std::endl;
    
    // for reading object
    const std::string readParamPath = "pdfparam.mpac";
    TlMsgPack readMsgPack;
    if (isVerbose == true) {
        std::cout << "read parmeter file: " << readParamPath << std::endl;
    }
    if (TlFile::isExist(readParamPath) == false) {
        std::cerr << "file not found: " << readParamPath << std::endl;
        return 1;
    }
    readMsgPack.load(readParamPath);
    const TlSerializeData readData = readMsgPack.getSerializeData();
    const TlOrbitalInfo orbitalInfo(readData["coordinates"],
                                    readData["basis_sets"]);

    DfEriEngine engine;
    const int shellTypeP = orbitalInfo.getShellType(shellIndexP);
    const int shellTypeQ = orbitalInfo.getShellType(shellIndexQ);
    const int shellTypeR = orbitalInfo.getShellType(shellIndexR);
    const int shellTypeS = orbitalInfo.getShellType(shellIndexS);
    const DfEriEngine::Query queryPQ(0, 0, shellTypeP, shellTypeQ);
    const DfEriEngine::Query queryRS(0, 0, shellTypeR, shellTypeS);
    const DfEriEngine::CGTO_Pair PQ = engine.getCGTO_pair(orbitalInfo,
                                                          shellIndexP, shellIndexQ,
                                                          0.0);
    const DfEriEngine::CGTO_Pair RS = engine.getCGTO_pair(orbitalInfo,
                                                          shellIndexR, shellIndexS,
                                                          0.0);
    engine.calc(queryPQ, queryRS, PQ, RS);

    // const int maxStepsP = 2 * shellTypeP + 1;
    // const int maxStepsQ = 2 * shellTypeQ + 1;
    // const int maxStepsR = 2 * shellTypeR + 1;
    // const int maxStepsS = 2 * shellTypeS + 1;
    const int maxStepsP = shellTypeP * (shellTypeP + 3) / 2 + 1;
    const int maxStepsQ = shellTypeQ * (shellTypeQ + 3) / 2 + 1;
    const int maxStepsR = shellTypeR * (shellTypeR + 3) / 2 + 1;
    const int maxStepsS = shellTypeS * (shellTypeS + 3) / 2 + 1;
    int index = 0;
    for (int i = 0; i < maxStepsP; ++i) {
        for (int j = 0; j < maxStepsQ; ++j) {
            for (int k = 0; k < maxStepsR; ++k) {
                for (int l = 0; l < maxStepsS; ++l) {
                    const double value = engine.WORK[index];
                    std::cout << TlUtils::format("(%2d %2d|%2d %2d)=%18.10f",
                                                 i, j, k, l, value)
                              << std::endl;
                    ++index;
                }
            }
        }
    }
    
    return EXIT_SUCCESS;
}


