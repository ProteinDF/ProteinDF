#include <iostream>
#include <cstdlib>
#include <vector>

#include "DfObject.h"
#include "TlGetopt.h"
#include "TlUtils.h"
#include "TlFile.h"
#include "TlOrbitalInfo.h"
#include "TlSymmetricMatrix.h"
#include "TlMsgPack.h"
#include "TlSerializeData.h"

#define AU_BOHR 1.889762

struct LCAO {
public:
    int index;
    double value;
};

struct LCAO_cmp {
public:
    bool operator()(const LCAO& a, const LCAO& b) {
        return (std::fabs(a.value) > std::fabs(b.value));
    }
};


void showHelp()
{
    std::cout << "pdftop <OPTIONS> MO_LEVEL NUM_OF_ITEMS" << std::endl;
    std::cout << " or " << std::endl;
    std::cout << "pdftop -A ATOM_INDEX" << std::endl;
    std::cout << std::endl;
    std::cout << "OPTIONS:" << std::endl;
    std::cout << "  <MO mode>" << std::endl;
    std::cout << "  -a:             output all basis functions" << std::endl;
    std::cout << "  -d DIRECTORY    set directory to read" << std::endl;
    std::cout << "  -h:             show help(this)" << std::endl;
    std::cout << "  -v:             enable verbose mode" << std::endl;
    std::cout << "  <atom mode>" << std::endl;
    std::cout << "  -A ATOM_INDEX   output MO levels in which AOs are major" << std::endl;
    std::cout << "  -T THRESHOLD    threshold of w_atom" << std::endl;
}


int main(int argc, char* argv[])
{
    TlGetopt opt(argc, argv, "ad:hr:vxA:T:");

    // parameters
    const bool isVerbose = (opt["v"] == "defined");
    std::string readDir = "./";
    if (opt["d"].empty() == false) {
        readDir = opt["d"];
    }

    bool isAllAOs = false;
    if (opt["a"] == "defined") {
        isAllAOs = true;
    }
    
    int level = 0;
    if (opt[1].empty() == false) {
        level = std::atoi(opt[1].c_str());
    }

    int top = 10;
    if (opt[2].empty() == false) {
        top = std::atoi(opt[2].c_str());
    }
    
    if (opt["h"] == "defined") {
        showHelp();
        return EXIT_FAILURE;
    }
    
    if (isVerbose == true) {
        std::cerr << "read dir: " << readDir << std::endl;
        std::cerr << "MO level: " << level << std::endl;
    }

    bool is_MO_mode = true;
    int inputAtomIndex = 0;
    if (opt["A"].empty() == false) {
        inputAtomIndex = std::atoi(opt["A"].c_str());
        is_MO_mode = false;
    }
    double threshold = 0.5;
    if (opt["T"].empty() == false) {
        threshold = std::atof(opt["T"].c_str());
    }

    // for reading object
    const std::string readParamPath = TlUtils::format("%s/pdfparam.mpac", readDir.c_str());
    TlMsgPack readMsgPack;
    readMsgPack.load(readParamPath);
    TlSerializeData pdfparam = readMsgPack.getSerializeData();
    const TlOrbitalInfo readOrbInfo(pdfparam["coordinates"],
                                    pdfparam["basis_sets"]);
    const int lastIteration = pdfparam["num_of_iterations"].getInt();
    DfObject dfObj(&pdfparam);

    // load S
    TlSymmetricMatrix S;
    S.load(dfObj.getSpqMatrixPath());
    
    // load C
    TlMatrix C;
    const std::string readMatrixPath = dfObj.getCMatrixPath(DfObject::RUN_RKS, lastIteration);
    if (isVerbose == true) {
        std::cout << TlUtils::format("read %s.", readMatrixPath.c_str()) << std::endl;
    }
    bool isRead = C.load(readMatrixPath);
    if (isVerbose == true) {
        std::cout << TlUtils::format("row = %d, col = %d", C.getNumOfRows(), C.getNumOfCols()) << std::endl;
    }
    if (C.getNumOfCols() <= level) {
        std::cerr << TlUtils::format("input level is out of range. input=%d, max=%d", level, C.getNumOfCols())
                  << std::endl;
        return EXIT_FAILURE;
    }
    const int numOfAOs = C.getNumOfRows();
    const int numOfMOs = C.getNumOfCols();

    const TlMatrix SC = S * C;
    if (is_MO_mode) {
        // MO mode -------------------------------------------------------------
        // double energy = 0.0;
        // {
        //     TlVector eigval;
        //     const std::string eigvalPath = dfObj.getEigenvaluesPath(DfObject::RUN_RKS, lastIteration);
        //     eigval.load(eigvalPath);
        //     energy = eigval.get(level);
        // }
        
        // MO
        const TlVector MO = C.getColVector(level);
        std::vector<LCAO> AOs(numOfAOs);
        double sum = 0.0;
        for (int i = 0; i < numOfAOs; ++i) {
            AOs[i].index = i;
            const double w = std::fabs(MO[i] * SC.get(i, level));
            sum += w;
            AOs[i].value = w;
        }
        
        // sort
        std::sort(AOs.begin(), AOs.end(), LCAO_cmp());
        
        // display
        const char shellTypes[] = {'s', 'p', 'd'};
        const char* angularMomentumTypes[] = {
            "s",
            "px", "py", "pz",
            "dxy", "dyz", "dzx", "dxx-yy", "dzz"
        };
        
        std::cout << TlUtils::format("MO level: %d", level)
                  << std::endl;
        // std::cout << TlUtils::format("energy = % 18.8f", energy)
        //           << std::endl;
        
        if (isAllAOs == true) {
            top = numOfAOs;
        }
        
        double accumulate = 0.0;
        for (int rank = 0; rank < top; ++rank) {
            const int AO = AOs[rank].index;
            const double w = AOs[rank].value;
            const TlPosition pos = readOrbInfo.getPosition(AO);
            
            const double rate = w / sum * 100.0;
            accumulate += rate;
            std::cout << TlUtils::format("No.%2d %5d %2s (% 8.3f, % 8.3f, % 8.3f) %c(%6s) coef=% 8.3f w=% 8.3f(%6.2f%%; %6.2f%%)",
                                         rank+1,
                                         readOrbInfo.getAtomIndex(AO),
                                         readOrbInfo.getAtomName(AO).c_str(),
                                         pos.x() / AU_BOHR, pos.y() / AU_BOHR, pos.z() / AU_BOHR,
                                         shellTypes[readOrbInfo.getShellType(AO)],
                                         angularMomentumTypes[readOrbInfo.getBasisType(AO)],
                                         MO.get(AO), w,
                                         rate, accumulate)
                      << std::endl;
        }
    } else {
        // 
        if (isVerbose) {
            std::cerr << TlUtils::format("threshold: %8.2e", threshold) << std::endl;
        }

        std::vector<int> pickupMO;
        for (int MO = 0; MO < numOfMOs; ++MO) {
            const TlVector MO_vec = C.getColVector(MO);

            if (isVerbose) {
                std::cerr << TlUtils::format("check %5dth MO ...", MO) << std::endl;
            }
            double sum = 0.0;
            double w_atom = 0.0;
            for (int i = 0; i < numOfAOs; ++i) {
                const double w = std::fabs(MO_vec[i] * SC.get(i, MO));
                sum += w;

                const int atomIndex = readOrbInfo.getAtomIndex(i);
                if (atomIndex == inputAtomIndex) {
                    w_atom += w;
                }
            }

            const double ratio = w_atom / sum;
            if (ratio > threshold) {
                std::cout << TlUtils::format("MO=%5d w_atom=% 8.3f", 
                                             MO, w_atom)
                          << std::endl;
            }
        }
    }

    return EXIT_SUCCESS;
}


