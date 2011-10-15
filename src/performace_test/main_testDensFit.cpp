#include <iostream>
#include <string>

#include "DfDensityFitting.h"

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
    log.setFilePath("ptDensFit.out");
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
    TlMsgPack mpac;
    mpac.load(mpacPath);
    TlSerializeData param = mpac.getSerializeData();
    TlOrbitalInfo orbitalInfo(param["coordinates"],
                              param["basis_sets"]);
    const int numOfAOs = orbitalInfo.getNumOfOrbitals();
    if (isVerbose == true) {
        std::cerr << TlUtils::format("number of AOs: %d", numOfAOs)
                  << std::endl;
    }

    // calc
    loggerTime("density fitting start");
    DfDensityFittingX dfDensityFitting(&param);
    dfDensityFitting.exec();
    loggerTime("end");

    // save
    loggerTime("finish");
    
    return EXIT_SUCCESS;
}

