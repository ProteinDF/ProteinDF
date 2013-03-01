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


void showHelp()
{
    std::cout << "pdf-info-orb [OPTIONS] index ..." << std::endl;
    std::cout << std::endl;
    std::cout << " display orbital information." << std::endl;
    std::cout << " index:        Specify orbital indeces which begin from '0' to output." << std::endl;
    std::cout << "               Output all orbitals if no entry." << std::endl;
    std::cout << " -a:           Atom index mode: index means atom serial number." << std::endl;
    std::cout << " -p PDF_param: ProteinDF parameter file (default: pdfparam.mpac)" << std::endl;
    std::cout << " -w mpac_file: output MessagaPack file" << std::endl;
    std::cout << " -v:           verbose" << std::endl;
       
}


int main(int argc, char* argv[])
{
    TlGetopt opt(argc, argv, "ahp:vw:");
    const bool isAtomIndexMode = (opt["a"] == "defined");
    const bool isVerbose = (opt["v"] == "defined");
    const bool isShowHelp = (opt["h"] == "defined");
    std::string pdfparamPath = "pdfparam.mpac";
    if (opt["p"].empty() != true) {
        pdfparamPath = opt["p"];
    }
    std::string saveMpacPath = "";
    if (opt["w"].empty() != true) {
        saveMpacPath = opt["w"];
    }

    if (isShowHelp) {
        showHelp();
        std::exit(0);
    }
    
    TlMsgPack mpac;
    if (isVerbose) {
        std::cerr << TlUtils::format("reading %s ...", pdfparamPath.c_str());
    }
    mpac.load(pdfparamPath);
    TlSerializeData data = mpac.getSerializeData();
    TlOrbitalInfo orbInfo(data["coordinates"], data["basis_sets"]);

    const std::size_t numOfAOs = orbInfo.getNumOfOrbitals();
    if (isVerbose == true) {
        std::cerr << "total atoms: " << orbInfo.getNumOfAtoms() << std::endl;
        std::cerr << "total orbitals: " << numOfAOs << std::endl;
    }

    const int numOfArgs = opt.getCount() -1; // opt[0] means the name of this program.

    if (isAtomIndexMode == true) {
        if (numOfArgs == 0) {
            showHelp();
            std::exit(1);
        }

        for (int i = 0; i < numOfArgs; ++i) {
            const int atomIndex = std::atoi(opt[i +1].c_str());
            
            std::cout << TlUtils::format(">>>> atom index: %d", atomIndex) << std::endl;
            for (int aoIndex = 0; aoIndex < numOfAOs; ++aoIndex) {
                if (atomIndex == orbInfo.getAtomIndex(aoIndex)) {
                    printOrbInfo(orbInfo, aoIndex);
                }
            }
        }
    } else {
        std::vector<int> orbitals;
        if (numOfArgs == 0) {
            orbitals.resize(numOfAOs);
            for (int i = 0; i < numOfAOs; ++i) {
                orbitals[i] = i;
            }
        } else {
            for (int i = 0; i < numOfArgs; ++i) {
                const int orb = std::atoi(opt[i +1].c_str());
                orbitals.push_back(orb);
            }
        }
        const int numOfQueries = orbitals.size();

        if (saveMpacPath.empty()) {
            for (int i = 0; i < numOfQueries; ++i) {
                const int orb = orbitals[i];
                if (orb < numOfAOs) {
                    printOrbInfo(orbInfo, orb);
                }
            }
        } else {
            TlSerializeData output;
            for (int i = 0; i < numOfQueries; ++i) {
                const int orb = orbitals[i];
                if (orb < numOfAOs) {
                    TlSerializeData entry;

                    entry["atom_index"] = orbInfo.getAtomIndex(orb);
                    entry["atomic_number"] = TlAtom::getElementNumber(orbInfo.getAtomName(orb));
                    entry["symbol"] = orbInfo.getAtomName(orb);
                    const TlPosition pos = orbInfo.getPosition(orb);
                    TlSerializeData xyz;
                    xyz.pushBack(pos.x() / AU_BOHR);
                    xyz.pushBack(pos.y() / AU_BOHR);
                    xyz.pushBack(pos.z() / AU_BOHR);
                    entry["xyz"] = xyz;

                    entry["shell_type_id"] = orbInfo.getShellType(orb);
                    entry["shell_type"] = TlUtils::format("%c", shellTypes[orbInfo.getShellType(orb)]);
                    entry["basis_type_id"] = orbInfo.getBasisType(orb);
                    entry["basis_type"] = std::string(angularMomentumTypes[orbInfo.getBasisType(orb)]);
                    
                    const int contractions = orbInfo.getCgtoContraction(orb);
                    for (int c = 0; c < contractions; ++c) {
                        TlSerializeData pGTO;
                        pGTO["exp"] = orbInfo.getExponent(orb, c);
                        pGTO["coef"] = orbInfo.getCoefficient(orb, c);
                        entry["pGTOs"].pushBack(pGTO);
                    }

                    output[orb] = entry;
                }
            }
            
            TlMsgPack mpack(output);
            if (isVerbose) {
                std::cerr << TlUtils::format("output: %s", saveMpacPath.c_str()) << std::endl;
            }
            mpack.save(saveMpacPath);
        }
    }

    return EXIT_SUCCESS;
}


