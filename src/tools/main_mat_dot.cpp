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

void showHelp(const std::string& name)
{
    std::cout << TlUtils::format("%s [options] input_file_path1 input_file_path2 output_file_path", name.c_str()) << std::endl;
    std::cout << " OPTIONS:" << std::endl;
    std::cout << "  -h:      show help" << std::endl;
    std::cout << "  -v:      verbose" << std::endl;
}

int main(int argc, char* argv[])
{
    TlGetopt opt(argc, argv, "hv");
    
    if (opt["h"] == "defined") {
        showHelp(opt[0]);
        return EXIT_SUCCESS;
    }
    
    const bool bVerbose = (opt["v"] == "defined");

    if (opt.getCount() <= 1) {
        showHelp(opt[0]);
        return EXIT_FAILURE;
    }
    const std::string inputMatrixPath1 = opt[1];
    const std::string inputMatrixPath2 = opt[2];
    const std::string outputMatrixPath = opt[3];

    if (bVerbose == true) {
        std::cerr << "load matrix: " << inputMatrixPath1 << std::endl;
    }

    TlMatrix A;
    if (TlMatrix::isLoadable(inputMatrixPath1)) {
        A.load(inputMatrixPath1);
    } else if (TlSymmetricMatrix::isLoadable(inputMatrixPath1)) {
        TlSymmetricMatrix tmp;
        tmp.load(inputMatrixPath1);
        A = tmp;
    } else {
        std::cerr << "can not open file: " << inputMatrixPath1 << std::endl;
        return EXIT_FAILURE;
    }

    if (bVerbose == true) {
        std::cerr << "load matrix: " << inputMatrixPath2 << std::endl;
    }
    TlMatrix B;
    if (TlMatrix::isLoadable(inputMatrixPath2)) {
        B.load(inputMatrixPath2);
    } else if (TlSymmetricMatrix::isLoadable(inputMatrixPath2)) {
        TlSymmetricMatrix tmp;
        tmp.load(inputMatrixPath2);
        B = tmp;
    } else {
        std::cerr << "can not open file: " << inputMatrixPath2 << std::endl;
        return EXIT_FAILURE;
    }

    if (bVerbose == true) {
        std::cerr << "running..." << std::endl;
    }

    A.dot(B);
    
    if (bVerbose == true) {
        std::cerr << "save matrix: " << outputMatrixPath << std::endl;
    }
    if (outputMatrixPath != "") {
        A.save(outputMatrixPath);
    }

    return EXIT_SUCCESS;
}


