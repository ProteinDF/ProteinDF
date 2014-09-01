#include <iostream>
#include <string>

#include "DfEriX.h"

#include "Fl_Geometry.h"
#include "Fl_Gto.h"
#include "TlLogging.h"
#include "TlMsgPack.h"
#include "TlSerializeData.h"
#include "TlGetopt.h"
#include "TlOrbitalInfo.h"
#include "TlTime.h"

void logger(const std::string& str) 
{
    TlLogging& log = TlLogging::getInstance();
    log.info(str);
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
    TlLogging& log = TlLogging::getInstance();
    log.setFilePath("ptERI.out");
    loggerTime("start");
    
    const bool isVerbose = (opt["v"] == "defined");
    
    std::string mpacPath = "pdfparam.mpac";
    if (opt["p"].empty() == false) {
        mpacPath = opt["p"];
    }

    std::string densityMatrixPath = "P.mat";
    if (opt["d"].empty() == false) {
        densityMatrixPath = opt["d"];
    }

    std::string KMatrixPath = "K.mat";
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
    if (isVerbose == true) {
        std::cerr << TlUtils::format("num of AOs: %d", numOfAOs)
                  << std::endl;
    }

    if (isVerbose == true) {
        std::cerr << TlUtils::format("density matrix: %s", densityMatrixPath.c_str())
                  << std::endl;
    }
    TlSymmetricMatrix P;
    P.load(densityMatrixPath);
    if (isVerbose == true) {
        std::cerr << TlUtils::format("dimensions of density matrix: %d", P.getNumOfRows())
                  << std::endl;
    }
    assert(P.getNumOfRows() == numOfAOs);
    
    // calc
    loggerTime("ERI start");
    DfEriX dfEri(&param);
    TlSymmetricMatrix K(numOfAOs);
    dfEri.getK(P, &K);
    loggerTime("ERI end");

    // save
    K.save(KMatrixPath);
    loggerTime("finish");
    
    return EXIT_SUCCESS;
}




