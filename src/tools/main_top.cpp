#include <iostream>
#include <cstdlib>
#include <vector>

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
    std::cout << "pdfTop <OPTIONS> MO_LEVEL NUM_OF_ITEMS" << std::endl;
    std::cout << std::endl;
    std::cout << "OPTIONS:" << std::endl;
    std::cout << "  -a:             output all basis functions" << std::endl;
    std::cout << "  -d [DIRECTORY]: set directory to read" << std::endl;
    std::cout << "  -h:             show help(this)" << std::endl;
    std::cout << "  -v:             enable verbose mode" << std::endl;
}


int main(int argc, char* argv[])
{
    TlGetopt opt(argc, argv, "ad:hr:vx");

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

    // for reading object
    const std::string readParamPath = TlUtils::format("%s/pdfparam.mpac", readDir.c_str());
    TlMsgPack readMsgPack;
    readMsgPack.load(readParamPath);
    const TlSerializeData readData = readMsgPack.getSerializeData();
    const TlOrbitalInfo readOrbInfo(readData["coordinates"],
                                    readData["basis_set"]);
    const int lastIteration = readData["iterations"].getInt();

    // load S
    TlSymmetricMatrix S;
    S.load(readDir + "/fl_Work/fl_Mtr_Spq.matrix");
    
    // load C
    TlMatrix C;
    std::string readMatrixPath = readDir + "/fl_Work/fl_Mtr_C.matrix.rks" + TlUtils::xtos(lastIteration);
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

    // energy
    double energy = 0.0;
    {
        TlVector eigval;
        const std::string eigvalPath = readDir + "/fl_Work/fl_Vct_Eigval" + TlUtils::xtos(lastIteration);
        eigval.load(eigvalPath);
        energy = eigval.get(level);
    }

    // MO
    const TlVector MO = C.getColVector(level);
    const int maxAO = MO.getSize();
    std::vector<LCAO> AOs(maxAO);
    for (int i = 0; i < maxAO; ++i) {
        AOs[i].index = i;
        const double value = MO[i];
        AOs[i].value = value;
    }
    std::sort(AOs.begin(), AOs.end(), LCAO_cmp());
    if (isAllAOs == true) {
        top = maxAO;
    }

    // calc normalize factor
    double allSquareCoef = 0.0; // 係数の2乗値の合計
    for (int i = 0; i < maxAO; ++i) {
        const double coef_i = MO[i];
        for (int j = 0; j < i; ++j) {
            const double coef_j = MO[j];
            allSquareCoef += 2.0 * coef_i * coef_j * S.get(i, j);
        }
        allSquareCoef += coef_i * coef_i;
    }
    
    // display
    const char shellTypes[] = {'s', 'p', 'd'};
    const char* angularMomentumTypes[] = {
        "s",
        "px", "py", "pz",
        "dxy", "dyz", "dzx", "dxx-yy", "dzz"
    };

    std::cout << TlUtils::format("MO level: %d", level)
              << std::endl;
    std::cout << TlUtils::format("energy = % 18.8f", energy)
              << std::endl;

    double accumulateRatio = 0.0;
    for (int rank = 0; rank < top; ++rank) {
        const int AO = AOs[rank].index;
        const double value = AOs[rank].value;

        const TlPosition pos = readOrbInfo.getPosition(AO);
        double squareCoef = 0.0; // 係数の2乗値
        for (int i = 0; i < maxAO; ++i) {
            squareCoef += MO[i] * value * S.get(AO, i);
        }
        const double ratio = 100.0 * squareCoef / allSquareCoef;
        accumulateRatio += ratio;

        std::cout << TlUtils::format("No.%2d %5d %2s (% 8.3f, % 8.3f, % 8.3f) %c(%6s) value=% 8.3f(%5.2f%%; %5.2f%%)",
                                     rank+1,
                                     readOrbInfo.getAtomIndex(AO),
                                     readOrbInfo.getAtomName(AO).c_str(),
                                     pos.x() / AU_BOHR, pos.y() / AU_BOHR, pos.z() / AU_BOHR,
                                     shellTypes[readOrbInfo.getShellType(AO)],
                                     angularMomentumTypes[readOrbInfo.getBasisType(AO)],
                                     value,
                                     ratio, accumulateRatio)
                  << std::endl;
    }
    
    return EXIT_SUCCESS;
}


