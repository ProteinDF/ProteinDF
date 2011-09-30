#include <iostream>
#include <string>

#include "DfCalcGridX.h"
#include "DfFunctional_B3LYP.h"

#include "Fl_Geometry.h"
#include "Fl_Gto.h"

#include "TlMsgPack.h"
#include "TlSerializeData.h"
#include "TlGetopt.h"
#include "TlLogX.h"
#include "TlOrbitalInfo.h"
#include "TlTime.h"

void logger(const std::string& str) 
{
    TlLogX& log = TlLogX::getInstance();
    log << str;
}


void loggerTime(const std::string& str) 
{
    std::string out = str;
    int size = out.size();
    if (size > 0) {
        if (out[size -1] == '\n') {
            out.erase(size -1, 1);
        }

        const std::string timeStr = "[" + TlTime::getNow() + "]";
        TlUtils::pad(out, (72 - timeStr.length()), ' ');
        out += (timeStr + "\n");
        logger(out);
    }
}


int main(int argc, char *argv[])
{
    TlGetopt opt(argc, argv, "d:hp:vs:");

    // setup log
    TlLogX& log = TlLogX::getInstance();
    log.setFilePath("evalEri.out");
    loggerTime("start");
    
    
    const bool isVerbose = (opt["v"] == "defined");
    
    std::string mpacPath = "pdfparam.mpac";
    if (opt["p"].empty() == false) {
        mpacPath = opt["p"];
    }

    std::string densityMatrixPath = "P.mtx";
    if (opt["d"].empty() == false) {
        densityMatrixPath = opt["d"];
    }

    std::string KMatrixPath = "Kxc.mtx";
    if (opt["s"].empty() == false) {
        KMatrixPath = opt["s"];
    }
    
    if (isVerbose == true) {
        std::cerr << TlUtils::format("PDF parameter: %s", mpacPath.c_str())
                  << std::endl;
    }
    TlMsgPack mpac;
    mpac.load(mpacPath);
    TlSerializeData param = mpac.getSerializeData();
    TlOrbitalInfo orbitalInfo(param["coordinates"],
                              param["basis_set"]);
    const int numOfAOs = orbitalInfo.getNumOfOrbitals();
    std::cerr << TlUtils::format("number of AOs: %d", numOfAOs)
              << std::endl;

    if (isVerbose == true) {
        std::cerr << TlUtils::format("density matrix: %s", densityMatrixPath.c_str())
                  << std::endl;
    }
    TlSymmetricMatrix P;
    P.load(densityMatrixPath);
    std::cerr << TlUtils::format("dimensions of density matrix: %d", P.getNumOfRows())
              << std::endl;
    assert(P.getNumOfRows() == numOfAOs);
    
    // cakc
    loggerTime("calc grid start");
    DfCalcGridX dfCalcGrid(&param);
    DfFunctional_B3LYP b3lyp;
    TlSymmetricMatrix K(numOfAOs);
    dfCalcGrid.calcXCIntegForFockAndEnergy(P, &b3lyp, &K);
    loggerTime("calc grid end");

    // save
    K.save(KMatrixPath);
    loggerTime("finish");
    
    return EXIT_SUCCESS;
}




