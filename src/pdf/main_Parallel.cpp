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

#ifdef HAVE_CONFIG_H
#include "config.h"    // this file created by autotools
#endif // HAVE_CONFIG_H

#include <sys/types.h>
#include <signal.h>
#include <unistd.h>
#include <cstring>
#include <cstdlib>
#include <ios>
#include <iostream>
#include <fstream>

#include "ProteinDF_Parallel.h"
#include "TlCommunicate.h"
#include "TlUtils.h"
#include "TlGetopt.h"
//#include "TlTime.h"

#ifdef __FUJITSU
#define PDF_MAIN MAIN__
#else
#define PDF_MAIN main
#endif // __FUJITSU

#define HOSTNAME_LEN 256

// signal(INT)を受け取ったときにMPI_Abort()を発行するため
static void func_int(int signum)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int proc = rComm.getRank();
    std::cerr << TlUtils::format("[%d] signal: %d received.",
                                 proc, signum)
              << std::endl;
    rComm.abort(1);
    std::abort();
}

// 例外処理時にMPI_Abort()を発行するため
void terminateHandler()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    rComm.abort(1);
    
    std::abort();
}

int PDF_MAIN(int argc, char *argv[])
{
    // for MPI
    TlCommunicate& rComm = TlCommunicate::getInstance(argc, argv);
    std::set_terminate(terminateHandler);
    
    // for signal
    struct sigaction act;
    memset(&act, 0, sizeof(act));
    act.sa_handler = func_int;
    if (sigaction(SIGINT, &act, NULL) < 0) {
        std::cerr << "signal handler could not registered." << std::endl;
    }
    
    // for debug
#ifndef NDEBUG
    {
        char* pHostName = new char[HOSTNAME_LEN];
        (void)gethostname(pHostName, HOSTNAME_LEN);
        const std::string hostname(pHostName);
        delete []pHostName;
        pHostName = NULL;

        const int numOfProcs = rComm.getNumOfProcs();
        for (int i = 0; i < numOfProcs; ++i) {
            if (i == rComm.getRank()) {
                std::cerr << TlUtils::format("PID %d on %s as rank %d",
                                             getpid(), hostname.c_str(), rComm.getRank())
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
#endif // NDEBUG

    // setup parameters
    TlGetopt opt(argc, argv, "Ddro:");

    bool isRestart = false;
    if (opt["r"] == "defined") {
        isRestart = true;
    }

    TlLogging& log = TlLogging::getInstance();
    std::string output = "fl_Out_Std";
    if (opt["o"].empty() != true) {
        output = opt["o"];
    }
    log.setFilePath(output);

    log.setProcID(rComm.getRank());
    if (opt["d"] == "defined") {
        log.setLevel(TlLogging::DEBUG, TlLogging::WARN);
    }
    if (opt["D"] == "defined") {
        log.setLevel(TlLogging::DEBUG, TlLogging::DEBUG);
    }
    
    // run ProteinDF
    ProteinDF_Parallel PDF;
    PDF.run();

    rComm.finalize();
    return EXIT_SUCCESS;
}

