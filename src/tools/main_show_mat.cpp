#include <iostream>
#include <cstdlib>

#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"
#include "TlGetopt.h"

int main(int argc, char* argv[])
{
    TlGetopt opt(argc, argv, "ghv");

    const bool bVerbose = (opt["v"] == "defined");
    const bool bGuessMode = (opt["g"] == "defined");

    std::string sPath = opt[1];
    if (bVerbose) {
        std::cerr << "loading... " << sPath << std::endl;
    }

    if (TlSymmetricMatrix::isLoadable(sPath) == true) {
        TlSymmetricMatrix M;
        M.load(sPath);

        if (bGuessMode == true) {
            M.saveText(std::cout);
        } else {
            M.print(std::cout);
        }
    } else if (TlMatrix::isLoadable(sPath) == true) {
        TlMatrix M;
        M.load(sPath);

        if (bGuessMode == true) {
            M.saveText(std::cout);
        } else {
            M.print(std::cout);
        }
    } else {
        std::cerr << "unknown file type: " << sPath << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}


