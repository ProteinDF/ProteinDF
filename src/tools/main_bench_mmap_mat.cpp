// Copyright (C) 2002-2014 The ProteinDF project
// see also AUTHORS and README.
//
// This file is part of ProteinDF.
//
// ProteinDF is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ProteinDF is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

#include <cstdlib>
#include <iostream>

#include <stdlib.h>
#include <time.h>

#include "TlMatrix.h"

#include "TlGetopt.h"
#include "TlMemManager.h"

int main(int argc, char* argv[]) {
    TlGetopt opt(argc, argv, "d:f:mrs");

    // std::size_t count = 10 * 1024 * 1024;

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
    TlDenseSymmetricMatrix_BLAS_Old S;
    S.load("S.matrix");
    TlDenseSymmetricMatrix_BLAS_Old eigVec;
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
