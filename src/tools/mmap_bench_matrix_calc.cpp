#include <iostream>
#include <cstdlib>

#include <stdlib.h>
#include <time.h>

#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"

#include "TlMemManager.h"
#include "TlGetopt.h"

int main(int argc, char* argv[])
{
    TlGetopt opt(argc, argv, "d:f:mrs");

    //std::size_t count = 10 * 1024 * 1024;

    std::string mmap_file = "/tmp/mmap_bench.map";
    if (!opt["f"].empty()) {
        mmap_file = opt["f"];
    }

    bool isMMap = false;
    if (!opt["m"].empty()) {
        isMMap = true;
    }

    // show information
    if (isMMap == true) {
        std::cout << "[mmap  ] ";
    } else {
        std::cout << "[memory] ";
    }

    // bench
    double prepared = 0.0;
    double elapse = 0.0;
    if (isMMap == true) {
        std::size_t dim = 15000;
        std::size_t needMemSize = dim * dim * sizeof(double) * 3;
        TlMemManager::setParam(needMemSize, mmap_file);
        TlMatrix::useMemManager(true);
    }

    clock_t startTime1 = clock();
    TlSymmetricMatrix S;
    S.load("S.matrix");
    TlSymmetricMatrix eigVec;
    TlVector eigVal;
    clock_t endTime1 = clock();
    
    clock_t startTime2 = clock();
    S.diagonal(&eigVal, &eigVec);
    clock_t endTime2 = clock();
    
    prepared = double(endTime1 - startTime1) / double(CLOCKS_PER_SEC);
    elapse = double(endTime2 - startTime2) / double(CLOCKS_PER_SEC);

    // results
    std::cout << "elapse time = " << elapse << " [s]" << std::endl;

    return EXIT_SUCCESS;
}


