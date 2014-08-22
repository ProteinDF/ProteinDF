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
#include "TlUtils.h"

void showHelp(const std::string& progname)
{
    std::cout << TlUtils::format("USAGE: %s [options] input_file_path output_file_path", progname.c_str()) << std::endl;
    std::cout << " OPTIONS:" << std::endl;
    std::cout << "  -h:      show help" << std::endl;
    std::cout << "  -v:      verbose" << std::endl;
}

int main(int argc, char* argv[])
{
    TlGetopt opt(argc, argv, "hvl:x:");
    
    if (opt["h"] == "defined") {
        showHelp(opt[0]);
        return EXIT_SUCCESS;
    }
    
    const bool bVerbose = (opt["v"] == "defined");

    if (opt.getCount() <= 1) {
        showHelp(opt[0]);
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


