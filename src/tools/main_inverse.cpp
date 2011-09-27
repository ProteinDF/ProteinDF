#include <iostream>
#include <cstdlib>

#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"
#include "TlGetopt.h"

void showHelp()
{
    std::cout << "inverse [options] input_file_path output_file_path" << std::endl;
    std::cout << " OPTIONS:" << std::endl;
    std::cout << "  -h:      show help" << std::endl;
    std::cout << "  -v:      verbose" << std::endl;
}

int main(int argc, char* argv[])
{
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
    const std::string inputMatrixPath = opt[1];
    const std::string outputMatrixPath = opt[2];

    if (bVerbose == true) {
        std::cerr << "load matrix: " << inputMatrixPath << std::endl;
    }
    if (TlMatrix::isLoadable(inputMatrixPath) == true) {
        TlMatrix A;
        A.load(inputMatrixPath);

        if (bVerbose == true) {
            std::cerr << "running..." << std::endl;
        }
        A.inverse();
    
        if (bVerbose == true) {
            std::cerr << "save inverse matrix: " << outputMatrixPath << std::endl;
        }
        if (outputMatrixPath != "") {
            A.save(outputMatrixPath);
        }
    } else if (TlMatrix::isLoadable(inputMatrixPath) != true) {
        TlSymmetricMatrix A;
        A.load(inputMatrixPath);

        if (bVerbose == true) {
            std::cerr << "running..." << std::endl;
        }
        A.inverse();
    
        if (bVerbose == true) {
            std::cerr << "save inverse matrix: " << outputMatrixPath << std::endl;
        }
        if (outputMatrixPath != "") {
            A.save(outputMatrixPath);
        }
    } else {
        std::cerr << "can not open file: " << inputMatrixPath << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}


