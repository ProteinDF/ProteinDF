#ifdef HAVE_CONFIG_H
#include "config.h"    // this file created by autotools
#endif // HAVE_CONFIG_H

#include <signal.h>
#include <cstring>
#include <cstdlib>
#include <ios>
#include <iostream>
#include <fstream>

#include "ProteinDF_Parallel.h"
#include "TlCommunicate.h"
#include "TlUtils.h"
#include "TlTime.h"

#ifdef AC_F77_MAIN
#define PDF_MAIN AC_F77_MAIN
#else
#define PDF_MAIN main
#endif // AC_F77_MAIN

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
    
    // run ProteinDF
    ProteinDF_Parallel PDF;
    PDF.run();

    rComm.finalize();
    return EXIT_SUCCESS;
}

