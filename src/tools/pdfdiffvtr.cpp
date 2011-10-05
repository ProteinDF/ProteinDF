#include <iostream>

#include "TlVector.h"
#include "TlGetopt.h"

void help();

int main(int argc, char* argv[])
{
    TlGetopt opt(argc, argv, "hv");

    bool bVerbose = (opt["v"] == "defined");
    if (opt.getCount() < 2) {
        help();
        std::exit(1);
    }

    const std::string sPath1 = opt[1];
    const std::string sPath2 = opt[2];
    if (bVerbose) {
        std::cerr << "loading... " << sPath1 << std::endl;
        std::cerr << "loading... " << sPath2 << std::endl;
    }

    int nErrorCode = 0;

    TlVector v1, v2;
    if (v1.load(sPath1) == false) {
        std::cerr << "could not open: " << sPath1 << std::endl;
        return 1;
    }
    if (v2.load(sPath2) == false) {
        std::cerr << "could not open: " << sPath2 << std::endl;
        return 1;
    }

    v1 -= v2;
    v1.print(std::cout);

    return nErrorCode;
}


void help()
{
    std::cout << "Usage: pdfdiffvtr [options]... FILE1 FILE2" << std::endl;
    std::cout << "compare ProteinDF vector files"
              << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << " -h           show help message (this)." << std::endl;
    std::cout << " -v           show message verbosely." << std::endl;
}


