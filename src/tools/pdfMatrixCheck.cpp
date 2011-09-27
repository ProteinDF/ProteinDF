#include <iostream>
#include <cstdlib>

#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"
#include "TlGetopt.h"

int main(int argc, char* argv[])
{
    TlGetopt opt(argc, argv, "ghv");

    const bool bVerbose = (opt["v"] == "defined");
    //const bool bGuessMode = (opt["g"] == "defined");

    std::string sPath = opt[1];
    if (bVerbose) {
        std::cerr << "loading... " << sPath << std::endl;
    }

    double threshold = 1.0;

    if (TlSymmetricMatrix::isLoadable(sPath) == true) {
        TlSymmetricMatrix M;
        M.load(sPath);

        const std::size_t maxIndex = M.getNumOfRows();
        for (std::size_t i = 0; i < maxIndex; ++i) {
            double value = std::fabs(M.get(i, i));

            if (value > threshold) {
                std::cout << TlUtils::format(" maxtrix(%6ld, %6ld) %e > (threshold) %e",
                                             i, i, value, threshold)
                          << std::endl;
            }
        }
    } else if (TlMatrix::isLoadable(sPath) == true) {
        TlMatrix M;
        M.load(sPath);

        const std::size_t maxIndex = std::min<std::size_t>(M.getNumOfRows(), M.getNumOfCols());
        for (std::size_t i = 0; i < maxIndex; ++i) {
            double value = std::fabs(M.get(i, i));

            if (value > threshold) {
                std::cout << TlUtils::format(" maxtrix(%6ld, %6ld) %e > (threshold) %e",
                                             i, i, value, threshold)
                          << std::endl;
            }
        }
    } else {
        std::cerr << "unknown file type: " << sPath << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}


