#include <iostream>

#include "TlVector.h"
#include "TlGetopt.h"

int main(int argc, char* argv[])
{
    TlGetopt opt(argc, argv, "hv");

    bool bVerbose = (opt["v"] == "defined");

    const std::string sPath1 = opt[1];
    const std::string sPath2 = opt[2];
    if (bVerbose) {
        std::cerr << "loading... " << sPath1 << std::endl;
        std::cerr << "loading... " << sPath2 << std::endl;
    }

    int nErrorCode = 0;

//   std::ifstream ifs1;
//   ifs1.open(sPath1.c_str());
//   if (ifs1.fail()){
//     std::cerr << "could not open file. " << sPath1 << std::endl;
//     return 1;
//   }

//   std::ifstream ifs2;
//   ifs2.open(sPath2.c_str());
//   if (ifs2.fail()){
//     std::cerr << "could not open file. " << sPath2 << std::endl;
//     return 1;
//   }

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


