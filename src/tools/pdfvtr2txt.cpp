#include <iostream>
#include <cstdlib>

#include "TlVector.h"
#include "TlGetopt.h"

void help();

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


void help()
{
    std::cout << "Usage: pdfvtr2txt [options]... FILE" << std::endl;
    std::cout << "display ProteinDF vector file"
              << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << " -g           show guess mode" << std::endl;
    std::cout << " -l           show list mode" << std::endl;
    std::cout << " -h           show help message (this)." << std::endl;
    std::cout << " -v           show message verbosely." << std::endl;
}


