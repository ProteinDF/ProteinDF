#include <iostream>
#include <cstdlib>
#include <vector>

#include "TlGetopt.h"
#include "TlOrbitalInfo.h"
#include "TlMsgPack.h"
#include "TlSerializeData.h"
#include "TlPseudoYaml.h"

#define AU_BOHR 1.889762

const char shellTypes[] = {'s', 'p', 'd'};
const char* angularMomentumTypes[] = {
    "s",
    "px", "py", "pz",
    "dxy", "dyz", "dzx", "dxx-yy", "dzz"
};


void printOrbInfo(const TlOrbitalInfo& tlOrbInfo,
                  std::size_t orb)
{
    const TlPosition pos = tlOrbInfo.getPosition(orb);
    std::cout << TlUtils::format("%d th orbital information:", orb)
              << std::endl;
    std::cout << TlUtils::format(" atom: index=%d, symbol=%s, position(% 8.3f, % 8.3f, % 8.3f)",
                                 tlOrbInfo.getAtomIndex(orb),
                                 tlOrbInfo.getAtomName(orb).c_str(),
                                 pos.x() / AU_BOHR, pos.y() / AU_BOHR, pos.z() / AU_BOHR)
              << std::endl;
    std::cout << TlUtils::format(" orbital: type=%c(%s)",
                                 shellTypes[tlOrbInfo.getShellType(orb)],
                                 angularMomentumTypes[tlOrbInfo.getBasisType(orb)])
              << std::endl;
    const int contractions = tlOrbInfo.getCgtoContraction(orb);
    std::cout << TlUtils::format("          contraction=%d",
                                 contractions)
              << std::endl;
    for (int j = 0; j < contractions; ++j) {
        std::cout << TlUtils::format("          coef=%e, exp=%e",
                                     tlOrbInfo.getCoefficient(orb, j),
                                     tlOrbInfo.getExponent(orb, j))
                  << std::endl;
    }
}


int main(int argc, char* argv[])
{
    TlGetopt opt(argc, argv, "ad:mv");
    const bool isAtomIndexMode = (opt["a"] == "defined");
    const bool isSaveMsgPack = (opt["m"] == "defined");
    const bool isVerbose = (opt["v"] == "defined");

    std::string readDir = ".";
    if (opt["d"].empty() == false) {
        readDir = opt["d"];
    }
    if (isVerbose == true) {
        std::cerr << "directory: " << readDir << std::endl;
    }

    TlMsgPack mpac;
    mpac.load("pdfparam.mpac");
    TlSerializeData data = mpac.getSerializeData();
    TlOrbitalInfo tlOrbInfo(data["model"]["coordinates"], data["model"]["basis_set"]);

    const std::size_t numOfAOs = tlOrbInfo.getNumOfOrbitals();
    if (isVerbose == true) {
        std::cerr << "total atoms: " << tlOrbInfo.getNumOfAtoms() << std::endl;
        std::cerr << "total orbitals: " << numOfAOs << std::endl;
    }

    const int numOfArgs = opt.getCount() -1; // opt[0] means the name of this program.
    if (isAtomIndexMode == true) {
        for (int i = 0; i < numOfArgs; ++i) {
            const int atomIndex = std::atol(opt[i +1].c_str());
            
            std::cout << TlUtils::format(">>>> atom index: %d", atomIndex) << std::endl;
            for (std::size_t aoIndex = 0; aoIndex < numOfAOs; ++aoIndex) {
                if (atomIndex == tlOrbInfo.getAtomIndex(aoIndex)) {
                    printOrbInfo(tlOrbInfo, aoIndex);
                }
            }
        }
    } else {
        if (isSaveMsgPack == true) {
            TlSerializeData output;
            for (int i = 0; i < numOfArgs; ++i) {
                const std::size_t orb = std::atol(opt[i +1].c_str());
                if (orb < numOfAOs) {
                    TlSerializeData tmp;
                    
                    tmp["atomic_number"] = TlAtom::getElementNumber(tlOrbInfo.getAtomName(orb));
                    const TlPosition pos = tlOrbInfo.getPosition(orb);
                    TlSerializeData coord;
                    coord.pushBack(pos.x());
                    coord.pushBack(pos.y());
                    coord.pushBack(pos.z());
                    tmp["coord"] = coord;
                    
                    output[orb] = tmp;
                }
            }
            
            TlMsgPack mpack(output);
            mpack.save("orbinfo.mpac");
        } else {
            for (int i = 0; i < numOfArgs; ++i) {
                const std::size_t orb = std::atol(opt[i +1].c_str());
                if (orb < numOfAOs) {
                    printOrbInfo(tlOrbInfo, orb);
                }
            }
        }
    }

    return EXIT_SUCCESS;
}


