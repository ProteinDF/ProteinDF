#include <mpi.h>

#ifdef _OPENMP
#include <omp.h>
#endif  // _OPENMMP

#include <cstdlib>
#include <iostream>
#include <string>

#include "TlCommunicate.h"
#include "TlUtils.h"

int main(int argc, char** argv) {
    TlCommunicate& rComm = TlCommunicate::getInstance(argc, argv);

    const int numOfProcs = rComm.getNumOfProcs();
    if (rComm.isMaster()) {
        std::cout << TlUtils::format("# procs: %d", numOfProcs) << std::endl;
    }
    rComm.barrier();

    for (int rank = 0; rank < numOfProcs; ++rank) {
        if (rank == rComm.getRank()) {
            const std::string hostName = rComm.getHostName();
            std::cout << TlUtils::format("[%d] %s", rank, hostName.c_str()) << std::endl;
        }
        rComm.barrier();
    }

    if (rComm.isMaster()) {
        std::string info = "";
#ifdef _OPENMP
        {
            const std::string ompInfo = TlUtils::format(" OpenMP threads: %d\n", omp_get_max_threads());
            info += ompInfo;

            {
                std::string ompProcBindStr = " OMP_PROC_BIND: ";
                const omp_proc_bind_t ompProcBind = omp_get_proc_bind();
                switch (ompProcBind) {
                    case omp_proc_bind_false:
                        ompProcBindStr += "False";
                        break;

                    case omp_proc_bind_true:
                        ompProcBindStr += "True";
                        break;

                    case omp_proc_bind_master:
                        ompProcBindStr += "Master";
                        break;

                    case omp_proc_bind_close:
                        ompProcBindStr += "Close";
                        break;

                    case omp_proc_bind_spread:
                        ompProcBindStr += "Spread";
                        break;

                    default:
                        ompProcBindStr += "UNKNOWN";
                }
                info += ompProcBindStr + "\n";
            }
        }
#else
        {
            info += " OpenMP is disabled.\n";
        }
#endif  // _OPENMP

        std::cout << info << std::endl;
    }

    rComm.finalize();
    return EXIT_SUCCESS;
}
