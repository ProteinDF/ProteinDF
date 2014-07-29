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

#include <iostream>
#include <cstdlib>

#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"
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

    std::string eigValPath = opt["l"];
    std::string eigVecPath = opt["x"];
    
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
    //const int numOfCols = A.getNumOfCols();

    TlMatrix eigVec(numOfDims, numOfDims);
    TlVector eigVal(numOfDims);
    
    if (bVerbose == true) {
        std::cerr << "running..." << inputMatrixPath << std::endl;
    }
    A.diagonal(&eigVal, &eigVec);
    
    if (bVerbose == true) {
        std::cerr << "save eigen values: " << eigValPath << std::endl;
    }
    if (eigValPath != "") {
        eigVal.save(eigValPath);
    }
    
    if (bVerbose == true) {
        std::cerr << "save eigen vectors: " << eigVecPath << std::endl;
    }
    if (eigVecPath != "") {
        eigVec.save(eigVecPath);
    }

    return EXIT_SUCCESS;
}


