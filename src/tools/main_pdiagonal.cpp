#include <iostream>
#include <cstdlib>

#include "TlCommunicate.h"
#include "TlSymmetricMatrix.h"
#include "TlDistributeMatrix.h"
#include "TlDistributeSymmetricMatrix.h"
#include "TlVector.h"
#include "TlGetopt.h"

void showHelp()
{
    std::cout << "diagonal [options] input_file_path" << std::endl;
    std::cout << " OPTIONS:" << std::endl;
    std::cout << "  -l FILE: save vector for eigen values" << std::endl;
    std::cout << "  -x FILE: save matrix for eigen vector" << std::endl;
    std::cout << "  -h:      show help" << std::endl;
    std::cout << "  -v:      verbose" << std::endl;
}

int main(int argc, char* argv[])
{
    TlCommunicate& rComm = TlCommunicate::getInstance(argc, argv);
    TlGetopt opt(argc, argv, "hvl:x:");
    
    if (opt["h"] == "defined") {
        showHelp();
        return EXIT_SUCCESS;
    }
    
    const bool bVerbose = (opt["v"] == "defined");

    if (opt.getCount() <= 1) {
        showHelp();
        return EXIT_FAILURE;
    }
    std::string inputMatrixPath = opt[1];

    std::string eigValPath = "";
    if (opt["l"].empty() == false) {
        eigValPath = opt["l"];
    }

    std::string eigVecPath = "";
    if (opt["x"].empty() == false) {
        eigVecPath = opt["x"];
    }
    
    if ((bVerbose == true) && (rComm.isMaster() == true)) {
        std::cerr << "load matrix: " << inputMatrixPath << std::endl;
    }
    if ((rComm.isMaster() == true) &&
        (TlSymmetricMatrix::isLoadable(inputMatrixPath) != true)) {
        std::cerr << "can not open file: " << inputMatrixPath << std::endl;
        return EXIT_FAILURE;
    }

    TlDistributeSymmetricMatrix A;
    A.load(inputMatrixPath);
    const int numOfDims = A.getNumOfRows();
    //const int numOfCols = A.getNumOfCols();

    TlDistributeMatrix eigVec(numOfDims, numOfDims);
    TlVector eigVal(numOfDims);
    
    if ((bVerbose == true) && (rComm.isMaster() == true)) {
        std::cerr << "running..." << inputMatrixPath << std::endl;
    }
    A.diagonal(&eigVal, &eigVec);
    
    if ((bVerbose == true) && (rComm.isMaster() == true)) {
        std::cerr << "save eigen values: " << eigValPath << std::endl;
    }
    if (eigValPath != "") {
        eigVal.save(eigValPath);
    }
    
    if ((bVerbose == true) && (rComm.isMaster() == true)) {
        std::cerr << "save eigen vectors: " << eigVecPath << std::endl;
    }
    if (eigVecPath != "") {
        eigVec.save(eigVecPath);
    }

    rComm.finalize();
    return EXIT_SUCCESS;
}


