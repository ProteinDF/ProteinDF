#include <iostream>
#include <cstdlib>

#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"
#include "TlGetopt.h"

void showHelp()
{
    std::cout << "cholesky [options] input_file_path" << std::endl;
    std::cout << " OPTIONS:" << std::endl;
    // std::cout << "  -l FILE: save vector for eigen values" << std::endl;
    // std::cout << "  -x FILE: save matrix for eigen vector" << std::endl;
    std::cout << "  -h:      show help" << std::endl;
    std::cout << "  -v:      verbose" << std::endl;
}

int main(int argc, char* argv[])
{
    TlGetopt opt(argc, argv, "hv");
    
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
    
    if (bVerbose == true) {
        std::cerr << "load matrix: " << inputMatrixPath << std::endl;
    }
    if (TlSymmetricMatrix::isLoadable(inputMatrixPath) != true) {
        std::cerr << "can not open file: " << inputMatrixPath << std::endl;
        return EXIT_FAILURE;
    }

    TlSymmetricMatrix A;
    A.load(inputMatrixPath);
    const int numOfDims = A.getNumOfRows();

    if (bVerbose == true) {
        std::cerr << "running..." << inputMatrixPath << std::endl;
    }
    TlMatrix L = A.choleskyFactorization2();

    TlMatrix Lt = L;
    Lt.transpose();

    TlMatrix LL = L * Lt;
    
    std::cout << ">>>> L" << std::endl;
    L.print(std::cout);
    std::cout << ">>>> A" << std::endl;
    A.print(std::cout);
    std::cout << ">>>> LL" << std::endl;
    LL.print(std::cout);
    
    return EXIT_SUCCESS;
}


