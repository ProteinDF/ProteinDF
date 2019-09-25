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

#include <cstdlib>
#include <iostream>

#include "DfRLMO.h"
#include "TlGetopt.h"
#include "TlMsgPack.h"
#include "TlOrbitalInfo.h"
#include "TlSerializeData.h"

void showHelp() {
    std::cout << "RLMO [options] start_block_AO_index ..." << std::endl;
    std::cout << " OPTIONS:" << std::endl;
    std::cout << "  -h:      show help" << std::endl;
    std::cout << "  -v:      verbose" << std::endl;
}

int main(int argc, char* argv[]) {
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
        int numOfBlocks = opt.getCount() - 1;
        for (int i = 0; i < numOfBlocks; ++i) {
            startBlockAOs.push_back(std::atoi(opt[i + 1].c_str()));
        }
    } else {
        TlOrbitalInfo orbInfo(param["coordinates"], param["basis_set"]);
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
