// Copyright (C) 2002-2014 The ProteinDF project
// see also AUTHORS and README.
//
// This file is part of ProteinDF.
//
// ProteinDF is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ProteinDF is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

#include <cstdlib>
#include <iostream>
#include <string>

#include "DfLocalize_Parallel.h"
#include "TlCommunicate.h"
#include "TlGetopt.h"
#include "TlLogging.h"
#include "TlMsgPack.h"

#define HOSTNAME_LEN 256

void showHelp(const std::string& progname) {
    std::cout << TlUtils::format("%s [OPTIONS] basisset_name ...", progname.c_str()) << std::endl;
    std::cout << std::endl;
    std::cout << " localize C matrix" << std::endl;
    std::cout << " -p PATH       set ProteinDF parameter file. default = "
                 "pdfparam.mpac"
              << std::endl;
    std::cout << " -c PATH       set C matrix path" << std::endl;
    std::cout << " -r            restart lo calculation" << std::endl;
    std::cout << " -h            show help" << std::endl;
    std::cout << " -v            verbose output" << std::endl;
}

int main(int argc, char* argv[]) {
    TlCommunicate& rComm = TlCommunicate::getInstance(argc, argv);

    TlGetopt opt(argc, argv, "c:dhp:rv");
    const bool isShowHelp = (opt["h"] == "defined");
    // const bool isVerbose = (opt["v"] == "defined");

// for debug
#ifndef NDEBUG
    {
        char* pHostName = new char[HOSTNAME_LEN];
        (void)gethostname(pHostName, HOSTNAME_LEN);
        const std::string hostname(pHostName);
        delete[] pHostName;
        pHostName = NULL;

        const int numOfProcs = rComm.getNumOfProcs();
        for (int i = 0; i < numOfProcs; ++i) {
            if (i == rComm.getRank()) {
                std::cerr << TlUtils::format("PID %d on %s as rank %d", getpid(), hostname.c_str(), rComm.getRank())
                          << std::endl;
            }
            rComm.barrier();
        }
        rComm.barrier();

        const int waitingTime = 5000;
        rComm.barrier();
        if (rComm.isMaster() == true) {
            std::cerr << TlUtils::format("waiting %d msec.", waitingTime) << std::endl;
        }
        TlTime::sleep(waitingTime);
        rComm.barrier();
    }
#endif  // NDEBUG

    if (isShowHelp) {
        if (rComm.isMaster() == true) {
            showHelp(opt[0]);
        }

        rComm.finalize();
        return EXIT_SUCCESS;
    }

    TlLogging& log = TlLogging::getInstance();
    log.setFilePath("lo-output.txt");

    if (opt["d"] == "defined") {
        log.setLevel(TlLogging::TL_DEBUG, TlLogging::TL_DEBUG);
    }

    std::string pdfParamPath = "pdfparam.mpac";
    if (opt["p"].empty() != true) {
        pdfParamPath = opt["p"];
    }

    std::string inputCMatrixPath = "";
    if (opt["c"].empty() != true) {
        inputCMatrixPath = opt["c"];
    }

    const bool isRestart = (opt["r"] == "defined");

    // setup
    TlSerializeData param;
    if (rComm.isMaster() == true) {
        if (opt["p"].empty() != true) {
            pdfParamPath = opt["p"];
        }

        TlMsgPack mpac;
        mpac.load(pdfParamPath);
        param = mpac.getSerializeData();
    }
    rComm.broadcast(param);

    // lo
    DfLocalize_Parallel lo(&param);
    if (!inputCMatrixPath.empty()) {
        lo.setCMatrixPath(inputCMatrixPath);
    }
    lo.setRestart(isRestart);

    lo.exec();

    if (rComm.isMaster() == true) {
        TlMsgPack mpac(param);
        mpac.save(pdfParamPath);
    }

    rComm.finalize();
    return EXIT_SUCCESS;
}
