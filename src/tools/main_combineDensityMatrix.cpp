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
    TlGetopt opt(argc, argv, "chr:vx");

    // parameters
    const bool isVerbose = (opt["v"] == "defined");
    const bool isCheck = (opt["c"] == "defined");
    const bool isCheckPosition = !(opt["x"] == "defined");
    double range = 1.0E-3;
    if (opt["r"].empty() == false) {
        range = std::atof(opt["r"].c_str());
    }
    
    std::string readDir = opt[1];
    std::string targetDir = opt[2];
    std::string savePpqPath = opt[3];
    if (savePpqPath == "") {
        savePpqPath = targetDir + "/fl_Work/fl_Mtr_Ppq.matrix.rks0";
    }

    if ((readDir == "") || (targetDir == "") ||
        (opt["h"] == "defined"))  {
        showHelp();
        return EXIT_FAILURE;
    }
    
    if (isVerbose == true) {
        std::cerr << "read dir: " << readDir << std::endl;
        std::cerr << "target dir: " << targetDir << std::endl;
        std::cerr << "save matrix path: " << savePpqPath << std::endl;
    }

    // for reading object
    const std::string readParamPath = TlUtils::format("%s/pdfparam.mpac", readDir.c_str());
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
    const TlOrbitalInfo readOrbInfo(readData["coordinates"],
                                    readData["basis_set"]);
    const int lastIteration = readData["iterations"].getInt();
    TlSymmetricMatrix readPpq;
    std::string readMatrixPath = readDir + "/fl_Work/fl_Mtr_Ppq.matrix.rks" + TlUtils::xtos(lastIteration);
    if (isVerbose == true) {
        std::cout << TlUtils::format("read %s.", readMatrixPath.c_str()) << std::endl;
    }
    bool isRead = readPpq.load(readMatrixPath);
    if (isRead == false) {
        std::cerr << "could not open: " << readMatrixPath << std::endl;
        return 2;
    }
    if (isVerbose == true) {
        std::cout << TlUtils::format("row = %d, col = %d", readPpq.getNumOfRows(), readPpq.getNumOfCols()) << std::endl;
    }
    
    // for Target
    const std::string targetParamPath = TlUtils::format("%s/pdfparam.mpac", targetDir.c_str());
    TlMsgPack targetMsgPack;
    targetMsgPack.load(targetParamPath);
    const TlSerializeData targetData = targetMsgPack.getSerializeData();
    const TlOrbitalInfo targetOrbInfo(targetData["coordinates"],
                                      targetData["basis_set"]);
    const int targetNumOfAOs = targetOrbInfo.getNumOfOrbitals();
    if (isVerbose == true) {
        std::cout << "target number Of AOs = " << targetNumOfAOs << std::endl;
    }
    TlSymmetricMatrix targetPpq(targetNumOfAOs);
    if (TlFile::isExist(savePpqPath) == true) {
        targetPpq.load(savePpqPath);
    }
    assert(targetPpq.getNumOfRows() == targetNumOfAOs);
    assert(targetPpq.getNumOfCols() == targetNumOfAOs);

    // Combine
    TlCombineDensityMatrix tlCombineDensMat(range, isVerbose);
    tlCombineDensMat.setCheckPosition(isCheckPosition);
    tlCombineDensMat.make(readOrbInfo, readPpq, targetOrbInfo, &targetPpq);

    targetPpq.save(savePpqPath);

    if (isCheck == true) {
        tlCombineDensMat.check(targetPpq);
    }

    return EXIT_SUCCESS;
}


