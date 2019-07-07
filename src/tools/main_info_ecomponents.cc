#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "TlGetopt.h"
#include "TlMsgPack.h"
#include "TlOrbitalInfo.h"
#include "TlSerializeData.h"
#include "TlUtils.h"

#include "df_total_energy_lapack.h"

int main(int argc, char* argv[]) {
    // parse args
    TlGetopt opt(argc, argv, "hi:m:p:v");
    const bool verbose = (opt["v"] == "defined");

    std::string pdfParamPath = "pdfparam.mpac";
    if (!opt["p"].empty()) {
        pdfParamPath = opt["p"];
    }

    // パラメータファイルの読み込み
    TlSerializeData param;
    {
        TlMsgPack mpac;
        mpac.load(pdfParamPath);
        param = mpac.getSerializeData();
    }

    //
    int iteration = -1;
    {
        if (!opt["i"].empty()) {
            iteration = std::atoi(opt["i"].c_str());
        } else {
            DfObject dfObj(&param);
            iteration = dfObj.iteration();
        }
    }
    if (verbose) {
        std::cerr << TlUtils::format("iteration: %d", iteration) << std::endl;
    }

    // prepare groups of atoms
    std::vector<std::vector<int> > groups_atom;

    // input by msgpack
    {
        TlSerializeData input;
        if (!opt["m"].empty()) {
            const std::string mpacPath = opt["m"];
            TlMsgPack mpac;
            mpac.load(mpacPath);
            input = mpac.getSerializeData();
        }

        TlSerializeData::ArrayIterator itEnd = input.endArray();
        for (TlSerializeData::ArrayIterator it = input.beginArray(); it != itEnd; ++it) {
            std::vector<int> atoms;
            TlSerializeData::ArrayIterator itAtomEnd = it->endArray();
            for (TlSerializeData::ArrayIterator itAtom = it->beginArray(); itAtom != itAtomEnd; ++itAtom) {
                atoms.push_back(itAtom->getInt());
            }
            groups_atom.push_back(atoms);
        }
    }

    // parse args
    {
        const int numOfArgs = opt.getCount();
        for (int i = 1; i < numOfArgs; ++i) {  // opt[0] is program name
            const std::string input = opt[i];
            std::cout << TlUtils::format("input %d: %s", i, input.c_str()) << std::endl;
            const std::vector<int> inputArray =
                TlUtils::nonreduntant_vector(TlUtils::vector_notation(input));
            groups_atom.push_back(inputArray);
        }
        if (groups_atom.size() == 0) {
            DfObject dfObj(&param);
            int numOfAtoms = dfObj.getNumOfAtoms();
            const std::string allAtoms = TlUtils::format("0-%d", numOfAtoms - 1);
            const std::vector<int> allAtomsArray =
              TlUtils::nonreduntant_vector(TlUtils::vector_notation(allAtoms));
            groups_atom.push_back(allAtomsArray);
        }
        if (verbose) {
            const int numOfGroups = groups_atom.size();
            for (int i = 0; i < numOfGroups; ++i) {
                std::cerr << TlUtils::format("group[%d] (atom): %s", i,
                                 TlUtils::vector2str(groups_atom[i]).c_str())
                          << std::endl;
            }
        }
    }

    DfTotalEnergy_Lapack dfTotalEnergy(&param);
    dfTotalEnergy.calc(iteration);
    dfTotalEnergy.output();

    //
    // double IE_total = 0.0;
    // const int numOfGroups = groups_atom.size();
    // for (int i = 0; i < numOfGroups; ++i) {
    //     {
    //         std::cout << TlUtils::format(">>>> group: %d <-> group: %d", i, i)
    //                   << std::endl;
    //         double IE = dfTotalEnergy.get_IE(groups_atom[i]);
    //         IE_total += IE;
    //     }
    //     for (int j = i + 1; j < numOfGroups; ++j) {
    //         std::cout << TlUtils::format(">>>> group: %d <-> group: %d", i, j)
    //                   << std::endl;
    //         // std::cout << TlUtils::format("group: %d [%s]", i,
    //         //                              TlUtils::vector2str(groups_atom[i]).c_str())
    //         //           << std::endl;
    //         // std::cout << TlUtils::format("group: %d [%s]", j,
    //         //                              TlUtils::vector2str(groups_atom[j]).c_str())
    //         //           << std::endl;
    //         double IE = dfTotalEnergy.get_IE(groups_atom[i], groups_atom[j]);
    //         IE_total += IE;
    //     }
    // }
    // std::cout << TlUtils::format("IE (total) = %lf", IE_total) << std::endl;
    dfTotalEnergy.get_IE(groups_atom);

    return EXIT_SUCCESS;
}
