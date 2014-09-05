#include <iostream>
#include <string>

#include "DfCalcGridX_Parallel.h"

#include "Fl_Geometry.h"
#include "Fl_Gto.h"

#include "TlMsgPack.h"
#include "TlSerializeData.h"
#include "TlGetopt.h"
#include "TlLogging.h"
#include "TlOrbitalInfo.h"
#include "TlTime.h"

#include "testXC_template.h"

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
    TlCommunicate& rComm = TlCommunicate::getInstance(argc, argv);
    TlGetopt opt(argc, argv, "d:hp:vs:");

    // setup log
    TlLogging& log = TlLogging::getInstance();
    log.setFilePath("testXC.out");
    loggerTime("start");
    
    const bool verbose = (opt["v"] == "defined");
    
    std::string mpacPath = "pdfparam.mpac";
    if (opt["p"].empty() == false) {
        mpacPath = opt["p"];
    }

    std::string densityMatrixPath = "P.mat";
    if (opt["d"].empty() == false) {
        densityMatrixPath = opt["d"];
    }

    std::string KMatrixPath = "Kxc.mat";
    if (opt["s"].empty() == false) {
        KMatrixPath = opt["s"];
    }
    
    if (verbose == true) {
        std::cerr << TlUtils::format("PDF parameter: %s", mpacPath.c_str())
                  << std::endl;
    }
    TlMsgPack mpac;
    mpac.load(mpacPath);
    TlSerializeData param = mpac.getSerializeData();

    // exec
    testXC<TlSymmetricMatrix, DfCalcGridX_Parallel>(param,
                                                    densityMatrixPath,
                                                    KMatrixPath,
                                                    verbose);

    rComm.finalize();
    return EXIT_SUCCESS;
}

