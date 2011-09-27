#include <iostream>
#include <cstdlib>

#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"
#include "TlGetopt.h"

void showHelp()
{
    std::cout << "multiple [options] input_file_path1 input_file_path2 output_file_path" << std::endl;
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
    const std::string inputMatrixPath1 = opt[1];
    const std::string inputMatrixPath2 = opt[2];
    const std::string outputMatrixPath = opt[3];

    if (bVerbose == true) {
        std::cerr << "load matrix: " << inputMatrixPath1 << std::endl;
    }
    if (TlMatrix::isLoadable(inputMatrixPath1) != true) {
        std::cerr << "can not open file: " << inputMatrixPath1 << std::endl;
        return EXIT_FAILURE;
    }

    TlMatrix A;
    A.load(inputMatrixPath1);

    if (bVerbose == true) {
        std::cerr << "load matrix: " << inputMatrixPath2 << std::endl;
    }
    if (TlMatrix::isLoadable(inputMatrixPath2) != true) {
        std::cerr << "can not open file: " << inputMatrixPath2 << std::endl;
        return EXIT_FAILURE;
    }

    TlMatrix B;
    B.load(inputMatrixPath2);

    if (bVerbose == true) {
        std::cerr << "running..." << std::endl;
    }

    TlMatrix C = A * B;
    
    if (bVerbose == true) {
        std::cerr << "save matrix: " << outputMatrixPath << std::endl;
    }
    if (outputMatrixPath != "") {
        C.save(outputMatrixPath);
    }

    return EXIT_SUCCESS;
}


