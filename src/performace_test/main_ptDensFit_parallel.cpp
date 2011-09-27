#include <iostream>
#include <string>

#include "DfDensityFittingX_Parallel.h"

#include "Fl_Geometry.h"
#include "Fl_Gto.h"

#include "TlCommunicate.h"
#include "TlMsgPack.h"
#include "TlSerializeData.h"
#include "TlGetopt.h"
#include "TlLogX.h"
#include "TlOrbitalInfo.h"
#include "TlTime.h"

void logger(const std::string& str) 
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        TlLogX& log = TlLogX::getInstance();
        log << str;
    }
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
    TlCommunicate& rComm = TlCommunicate::getInstance(argc, argv);
    TlGetopt opt(argc, argv, "d:hp:vs:");

    // setup log
    TlLogX& log = TlLogX::getInstance();
    log.setFilePath("evalDensFit.out");
    loggerTime("start");
    
    
    const bool isVerbose = (opt["v"] == "defined");
    
    std::string mpacPath = "pdfparam.mpac";
    if (opt["p"].empty() == false) {
        mpacPath = opt["p"];
    }

    if (isVerbose == true) {
        std::cerr << TlUtils::format("PDF parameter: %s", mpacPath.c_str())
                  << std::endl;
    }
    TlSerializeData param;
    if (rComm.isMaster() == true) {
        TlMsgPack mpac;
        mpac.load(mpacPath);
        param = mpac.getSerializeData();
    }
    rComm.broadcast(param);
    TlOrbitalInfo orbitalInfo(param["model"]["coordinates"],
                              param["model"]["basis_set"]);
    const int numOfAOs = orbitalInfo.getNumOfOrbitals();
    std::cerr << TlUtils::format("number of AOs: %d", numOfAOs)
              << std::endl;

    // cakc
    loggerTime("density fitting start");
    {
        DfDensityFittingX_Parallel dfDensityFitting(&param);
        dfDensityFitting.exec();
    }
    loggerTime("end");

    // save
    loggerTime("finish");

    std::cerr << "check..." << std::endl;
    rComm.checkNonBlockingCommunications();
    std::cerr << "done" << std::endl;
    
    rComm.finalize();
    return EXIT_SUCCESS;
}




