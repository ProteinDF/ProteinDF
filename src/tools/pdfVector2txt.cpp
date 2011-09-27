#include <iostream>
#include <cstdlib>

#include "TlVector.h"
#include "TlGetopt.h"

int main(int argc, char* argv[])
{
    TlGetopt opt(argc, argv, "ghlv");

    const bool bGuessMode = (opt["g"] == "defined");
    const bool bListMode = (opt["l"] == "defined");
    const bool bVerbose = (opt["v"] == "defined");

    std::string sPath = opt[1];
    if (bVerbose) {
        std::cerr << "loading... " << sPath << std::endl;
    }

    TlVector v;
    bool bIsLoadable = v.load(sPath);
    if (bIsLoadable == false) {
        std::cerr << "could not open: " << sPath << std::endl;
        return EXIT_FAILURE;
    }

    if (bGuessMode == true) {
        v.outputText(std::cout);
    } else if (bListMode == true) {
        const int nSize = v.getSize();
        for (int i = 0; i < nSize; ++i) {
            std::cout << v[i] << std::endl;
        }
    } else {
        v.print(std::cout);
    }

    return EXIT_SUCCESS;
}


