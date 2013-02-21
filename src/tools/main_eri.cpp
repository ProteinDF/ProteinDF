#include <iostream>
#include <cstdlib>
#include <vector>

#include "DfEriEngine.h"
#include "TlGetopt.h"
#include "TlUtils.h"
#include "TlOrbitalInfo.h"
#include "TlMsgPack.h"
#include "TlFile.h"

void showHelp()
{
    std::cout << "pdf-eri indexP indexQ indexR indexS" << std::endl;
    std::cout << "OPTIONS:" << std::endl;
    std::cout << "  -p PDF_PARAM: ProteinDF parameter file";
    std::cout << "  -h:           show help(this)" << std::endl;
    std::cout << "  -v:           verbose" << std::endl;
}


void calc_eri(const int shellIndexP, const int shellIndexQ,
              const int shellIndexR, const int shellIndexS,
              const TlOrbitalInfo& orbitalInfo);

void output(int indexP, int indexQ,
            int indexR, int indexS,
            const TlOrbitalInfoObject& orbitalInfo,
            const double* WORK);

int main(int argc, char* argv[])
{
    TlGetopt opt(argc, argv, "hp:v");

    // parameters
    const bool isVerbose = (opt["v"] == "defined");

    int indexP = std::atoi(opt[1].c_str());
    int indexQ = std::atoi(opt[2].c_str());
    int indexR = std::atoi(opt[3].c_str());
    int indexS = std::atoi(opt[4].c_str());
    
    // for reading object
    std::string readParamPath = "pdfparam.mpac";
    if (!opt["p"].empty()) {
        readParamPath = opt["p"];
    }
    TlMsgPack readMsgPack;
    if (isVerbose == true) {
        std::cout << "read parmeter file: " << readParamPath << std::endl;
    }
    if (TlFile::isExist(readParamPath) == false) {
        std::cerr << "file not found: " << readParamPath << std::endl;
        return 1;
    }
    readMsgPack.load(readParamPath);
    TlSerializeData readData = readMsgPack.getSerializeData();
    const TlOrbitalInfo orbitalInfo(readData["coordinates"],
                                    readData["basis_sets"]);

    const int numOfAOs = orbitalInfo.getNumOfOrbitals();
    if (! (0 <=indexP) && (indexP < numOfAOs)) {
        std::cerr << TlUtils::format("illegal input parameter1. %d / %d",
                                     indexP, numOfAOs) 
                  << std::endl;
        std::exit(1);
    }
    if (! (0 <=indexQ) && (indexQ < numOfAOs)) {
        std::cerr << TlUtils::format("illegal input parameter2. %d / %d",
                                     indexQ, numOfAOs) 
                  << std::endl;
        std::exit(1);
    }
    if (! (0 <=indexR) && (indexR < numOfAOs)) {
        std::cerr << TlUtils::format("illegal input parameter3. %d / %d",
                                     indexR, numOfAOs) 
                  << std::endl;
        std::exit(1);
    }
    if (! (0 <=indexS) && (indexS < numOfAOs)) {
        std::cerr << TlUtils::format("illegal input parameter4. %d / %d",
                                     indexS, numOfAOs) 
                  << std::endl;
        std::exit(1);
    }

    calc_eri(indexP, indexQ, indexR, indexS,
             orbitalInfo);

    // if (opt["o"] != "defined") {
    // } else {
    //     std::cerr << "old engine" << std::endl;
    //     calc_eri_old(&readData,
    //                  shellIndexP, shellIndexQ, shellIndexR, shellIndexS,
    //                  orbitalInfo);
    // }
    
    return EXIT_SUCCESS;
}


void calc_eri(const int indexP, const int indexQ,
              const int indexR, const int indexS,
              const TlOrbitalInfo& orbitalInfo)
{
    DfEriEngine engine;
    const int shellIndexP = orbitalInfo.getShellIndex(indexP);
    const int shellIndexQ = orbitalInfo.getShellIndex(indexQ);
    const int shellIndexR = orbitalInfo.getShellIndex(indexR);
    const int shellIndexS = orbitalInfo.getShellIndex(indexS);
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

    output(indexP, indexQ, indexR, indexS,
           orbitalInfo,
           engine.WORK);
}


void output(int indexP, int indexQ,
            int indexR, int indexS,
            const TlOrbitalInfoObject& orbitalInfo,
            const double* WORK)
{
    const int shellIndexP = orbitalInfo.getShellIndex(indexP);
    const int shellIndexQ = orbitalInfo.getShellIndex(indexQ);
    const int shellIndexR = orbitalInfo.getShellIndex(indexR);
    const int shellIndexS = orbitalInfo.getShellIndex(indexS);
    const int basisTypeP = indexP - shellIndexP;
    const int basisTypeQ = indexQ - shellIndexQ;
    const int basisTypeR = indexR - shellIndexR;
    const int basisTypeS = indexS - shellIndexS;
    const int shellTypeP = orbitalInfo.getShellType(shellIndexP);
    const int shellTypeQ = orbitalInfo.getShellType(shellIndexQ);
    const int shellTypeR = orbitalInfo.getShellType(shellIndexR);
    const int shellTypeS = orbitalInfo.getShellType(shellIndexS);
    const int maxStepsP = 2 * shellTypeP + 1;
    const int maxStepsQ = 2 * shellTypeQ + 1;
    const int maxStepsR = 2 * shellTypeR + 1;
    const int maxStepsS = 2 * shellTypeS + 1;
    
    {
        int index = 0;
        for (int i = 0; i < maxStepsP; ++i) {
            const int P = shellIndexP + i;
            for (int j = 0; j < maxStepsQ; ++j) {
                const int Q = shellIndexQ + j;
                for (int k = 0; k < maxStepsR; ++k) {
                    const int R = shellIndexR + k;
                    for (int l = 0; l < maxStepsS; ++l) {
                        const int S = shellIndexS + l;
                        const double value = WORK[index];
                        std::cout << TlUtils::format("[%d %d|%d %d] -> % f", P, Q, R, S, value)
                                  << std::endl;
                        ++index;
                    }
                }
            }
        }
    }

    const int index = ((basisTypeP * maxStepsQ + basisTypeQ) * maxStepsR + basisTypeR) * maxStepsS + basisTypeS;
    const double value = WORK[index];
    std::cout << TlUtils::format("(%2d %2d|%2d %2d)=%18.10f",
                                 indexP, indexQ, indexR, indexS, value)
              << std::endl;
}
