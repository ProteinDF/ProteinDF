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

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <string>

#include "DfPopulation.h"
#include "TlGetopt.h"
#include "TlMsgPack.h"
#include "TlSerializeData.h"
#include "tl_dense_vector_lapack.h"

void usage();

std::vector<std::vector<int> > getGroupsOfAtoms(const std::string& mpacPath) {
    // prepare groups of atoms
    std::vector<std::vector<int> > groups_atom;

    TlMsgPack mpac;
    mpac.load(mpacPath);
    TlSerializeData input = mpac.getSerializeData();

    TlSerializeData::ArrayIterator itEnd = input.endArray();
    for (TlSerializeData::ArrayIterator it = input.beginArray(); it != itEnd;
         ++it) {
        std::vector<int> atoms;
        TlSerializeData::ArrayIterator itAtomEnd = it->endArray();
        for (TlSerializeData::ArrayIterator itAtom = it->beginArray();
             itAtom != itAtomEnd; ++itAtom) {
            atoms.push_back(itAtom->getInt());
        }
        groups_atom.push_back(atoms);
    }

    return groups_atom;
}

int main(int argc, char* argv[]) {
    TlGetopt opt(argc, argv, "hi:m:vs:");

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
        iteration = param["num_of_iterations"].getInt();
    }

    std::vector<std::vector<int> > groups_atom;
    if (!opt["m"].empty()) {
        const std::string mpacPath = opt["m"];
        groups_atom = getGroupsOfAtoms(mpacPath);
    }

    if (isVerbose == true) {
        std::cerr << "iteration = " << iteration << std::endl;
    }

    DfPopulation dfPop(&param);
    dfPop.exec(iteration);
    dfPop.getReport(iteration, std::cout);

    // if (savePath.empty() != true) {
    //   // const TlDenseGeneralMatrix_Lapack mtx =
    //   dfPop.getAtomPopData(iteration); if (isVerbose == true) {
    //     std::cerr << "save Mulliken Data as " << savePath << std::endl;
    //   }
    //   //mtx.save(savePath);
    // } else {
    //   dfPop.getReport(iteration, std::cout);
    // }

    //
    int groupIndex = 0;
    std::vector<std::vector<int> >::const_iterator itGrpEnd = groups_atom.end();
    for (std::vector<std::vector<int> >::const_iterator itGrp =
             groups_atom.begin();
         itGrp != itGrpEnd; ++itGrp) {
        double grpSum = 0.0;
        std::vector<int>::const_iterator itEnd = itGrp->end();
        for (std::vector<int>::const_iterator it = itGrp->begin(); it != itEnd;
             ++it) {
            const double charge = dfPop.getCharge(*it);
            grpSum += charge;
        }
        std::cout << TlUtils::format("group [%4d]: % 8.3f", groupIndex, grpSum)
                  << std::endl;
        ++groupIndex;
    }

    return EXIT_SUCCESS;
}

void usage() {
    std::cerr << "usage: pdf-pop-mulliken [options] " << std::endl;
    std::cerr << std::endl;
    std::cerr << "calculation Mulliken Population." << std::endl;
    std::cerr << "  -p param:\n";
    std::cerr << "  -i num:    set SCF iteration to get Mulliken Population\n";
    std::cerr << "  -s path:   save atom population data as pdf matrix file\n";
    std::cerr << "  -h:        show help(this)\n";
    std::cerr << "  -v:        verbose\n";
}
