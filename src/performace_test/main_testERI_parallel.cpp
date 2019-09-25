#include <cstdlib>
#include <iostream>
#include <string>

#include "DfEriX_Parallel.h"
#include "TlCommunicate.h"
#include "TlGetopt.h"
#include "TlLogging.h"
#include "TlMsgPack.h"
#include "TlOrbitalInfo.h"
#include "TlSerializeData.h"

int main(int argc, char* argv[]) {
    TlCommunicate& rComm = TlCommunicate::getInstance(argc, argv);

    TlGetopt opt(argc, argv, "d:hp:vs:");

    const bool isVerbose = (opt["v"] == "defined");

    std::string mpacPath = "pdfparam.mpac";
    if (opt["p"].empty() == false) {
        mpacPath = opt["p"];
    }

    std::string densityMatrixPath = "P.mtx";
    if (opt["d"].empty() == false) {
        mpacPath = opt["d"];
    }

    std::string KMatrixPath = "K.mtx";
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
    TlOrbitalInfo orbitalInfo(param["coordinates"], param["basis_set"]);
    const int numOfAOs = orbitalInfo.getNumOfOrbitals();
    if (isVerbose == true) {
        std::cerr << TlUtils::format("num of AOs: %d", numOfAOs) << std::endl;
    }

    if (isVerbose == true) {
        std::cerr << TlUtils::format("density matrix: %s",
                                     densityMatrixPath.c_str())
                  << std::endl;
    }
    TlDistributeSymmetricMatrix P;
    P.load(densityMatrixPath);
    if (isVerbose == true) {
        std::cerr << TlUtils::format("dimensions of density matrix: %d",
                                     P.getNumOfRows())
                  << std::endl;
    }
    assert(P.getNumOfRows() == numOfAOs);

    // setup log
    TlLogging& log = TlLogging::getInstance();
    log.setFilePath("evalEri.out");

    // cakc
    DfEriX_Parallel dfEri(&param);
    TlDistributeSymmetricMatrix K(numOfAOs);
    dfEri.getK_D(P, &K);

    // save
    K.save(KMatrixPath);

    rComm.finalize();
    return EXIT_SUCCESS;
}
