#include <iostream>
#include <cstdlib>

#include "TlMatrix.h"
#include "TlGetopt.h"

void showHelp()
{
    std::cout << "transpose [options] input_path output_path" << std::endl;
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

    if (opt.getCount() <= 2) {
        showHelp();
        return EXIT_FAILURE;
    }
    std::string inputMatrixPath = opt[1];
    std::string outputMatrixPath = opt[2];

    if (bVerbose == true) {
        std::cerr << "load matrix: " << inputMatrixPath << std::endl;
    }

    TlMatrix A;
    A.load(inputMatrixPath);

    A.transpose();
    
    if (bVerbose == true) {
        std::cerr << "save matrix: " << outputMatrixPath << std::endl;
    }
    A.save(outputMatrixPath);
    
    return EXIT_SUCCESS;
}


