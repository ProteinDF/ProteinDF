#include <iostream>
#include <cstdlib>

#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"
#include "TlGetopt.h"

int main(int argc, char* argv[])
{
    TlGetopt opt(argc, argv, "ghv");

//   const bool bVerbose = (opt["v"] == "defined");
//   const bool bGuessMode = (opt["g"] == "defined");

    std::string sPath1 = opt[1];
    std::string sPath2 = opt[2];
    std::string sPath3 = opt[3];
//   if (bVerbose){
//     std::cerr << "loading... " << sPath << std::endl;
//     std::cerr << "loading... " << sPath << std::endl;
//   }

    if (TlSymmetricMatrix::isLoadable(sPath1) != true) {
        std::cerr << "can not open file: " << sPath1 << std::endl;
        return EXIT_FAILURE;
    }
    if (TlSymmetricMatrix::isLoadable(sPath2) != true) {
        std::cerr << "can not open file: " << sPath1 << std::endl;
        return EXIT_FAILURE;
    }

    TlSymmetricMatrix M1;
    M1.load(sPath1);
    //const int nRows1 = M1.getNumOfRows();
    //const int nCols1 = M1.getNumOfCols();

    TlSymmetricMatrix M2;
    M2.load(sPath2);
    //const int nRows2 = M2.getNumOfRows();
    //const int nCols2 = M2.getNumOfCols();

    M1 += M2;

    M1.save(sPath3);

    return EXIT_SUCCESS;
}


