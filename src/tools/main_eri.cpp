// Copyright (C) 2002-2014 The ProteinDF project
// see also AUTHORS and README.
// 
// This file is part of ProteinDF.
// 
// ProteinDF is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// ProteinDF is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

#include <iostream>
#include <cstdlib>
#include <vector>

#include "DfEriEngine.h"
#include "TlGetopt.h"
#include "TlUtils.h"
#include "TlOrbitalInfo.h"
#include "TlMsgPack.h"
#include "TlFile.h"

struct ShellIndex4 {
    int p;
    int q;
    int r;
    int s;
};

std::vector<ShellIndex4> getList(const std::string& listPath);

void calc_eri(DfEriEngine* pEngine,
              const int shellIndexP, const int shellIndexQ,
              const int shellIndexR, const int shellIndexS,
              const TlOrbitalInfo& orbitalInfo);

void output(int indexP, int indexQ,
            int indexR, int indexS,
            const TlOrbitalInfoObject& orbitalInfo,
            const double* WORK);

void showHelp();

int main(int argc, char* argv[])
{
    TlGetopt opt(argc, argv, "hp:v");

    // parameters
    const bool isVerbose = (opt["v"] == "defined");
    std::string readParamPath = "pdfparam.mpac";
    if (!opt["p"].empty()) {
        readParamPath = opt["p"];
    }

    if (opt["h"] == "defined") {
        showHelp();
        return EXIT_SUCCESS;
    }

    // for reading object
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
                                    readData["basis_set"]);
    const int numOfAOs = orbitalInfo.getNumOfOrbitals();

    // calc ERI
    if (opt.getCount() == 2) {
        std::string listPath = opt[1];
        std::vector<ShellIndex4> indeces = getList(listPath);
        const int size = indeces.size();
        DfEriEngine engine;
        for (int i = 0; i < size; ++i) {
            int indexP = indeces[i].p;
            int indexQ = indeces[i].q;
            int indexR = indeces[i].r;
            int indexS = indeces[i].s;
            
            calc_eri(&engine,
                     indexP, indexQ, indexR, indexS,
                     orbitalInfo);
        }
    } else if (opt.getCount() == 5) {
        int indexP = std::atoi(opt[1].c_str());
        int indexQ = std::atoi(opt[2].c_str());
        int indexR = std::atoi(opt[3].c_str());
        int indexS = std::atoi(opt[4].c_str());

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
        
        DfEriEngine engine;
        calc_eri(&engine,
                 indexP, indexQ, indexR, indexS,
                 orbitalInfo);
    } else {
        showHelp();
        return EXIT_SUCCESS;
    }
    

    return EXIT_SUCCESS;
}


std::vector<ShellIndex4> getList(const std::string& listPath)
{
    std::vector<ShellIndex4> answer;

    TlMsgPack mpac;
    mpac.load(listPath);
    TlSerializeData data = mpac.getSerializeData();
    if (data.getType() == TlSerializeData::ARRAY) {
        const int size = data.getSize();
        answer.resize(size);
        for (int i = 0; i < size; ++i) {
            const TlSerializeData& item = data.getAt(i);
            answer[i].p = item.getAt(0).getInt();
            answer[i].q = item.getAt(1).getInt();
            answer[i].r = item.getAt(2).getInt();
            answer[i].s = item.getAt(3).getInt();
        }
    }

    return answer;
}


void calc_eri(DfEriEngine* pEngine,
              const int indexP, const int indexQ,
              const int indexR, const int indexS,
              const TlOrbitalInfo& orbitalInfo)
{
    const int shellIndexP = orbitalInfo.getShellIndex(indexP);
    const int shellIndexQ = orbitalInfo.getShellIndex(indexQ);
    const int shellIndexR = orbitalInfo.getShellIndex(indexR);
    const int shellIndexS = orbitalInfo.getShellIndex(indexS);
    const int shellTypeP = orbitalInfo.getShellType(shellIndexP);
    const int shellTypeQ = orbitalInfo.getShellType(shellIndexQ);
    const int shellTypeR = orbitalInfo.getShellType(shellIndexR);
    const int shellTypeS = orbitalInfo.getShellType(shellIndexS);
    // const DfEriEngine::AngularMomentum2 queryPQ(0, 0, shellTypeP, shellTypeQ);
    // const DfEriEngine::AngularMomentum2 queryRS(0, 0, shellTypeR, shellTypeS);
    // const DfEriEngine::CGTO_Pair PQ = pEngine->getCGTO_pair(orbitalInfo,
    //                                                         shellIndexP, shellIndexQ,
    //                                                         0.0);
    // const DfEriEngine::CGTO_Pair RS = pEngine->getCGTO_pair(orbitalInfo,
    //                                                         shellIndexR, shellIndexS,
    //                                                         0.0);
    // pEngine->calc(queryPQ, queryRS, PQ, RS);
    pEngine->calc(0, orbitalInfo, shellIndexP,
                  0, orbitalInfo, shellIndexQ,
                  0, orbitalInfo, shellIndexR,
                  0, orbitalInfo, shellIndexS);

    output(indexP, indexQ, indexR, indexS,
           orbitalInfo,
           pEngine->WORK);
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

void showHelp()
{
    std::cout << "pdf-eri indexP indexQ indexR indexS" << std::endl;
    std::cout << " or " << std::endl;
    std::cout << "pdf-eri index_list" << std::endl;
    std::cout << "OPTIONS:" << std::endl;
    std::cout << "  -p PDF_PARAM: ProteinDF parameter file";
    std::cout << "  -h:           show help(this)" << std::endl;
    std::cout << "  -v:           verbose" << std::endl;
}


