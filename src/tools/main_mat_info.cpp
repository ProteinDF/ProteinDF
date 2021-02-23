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
#include <string>

#include "TlGetopt.h"
#include "TlUtils.h"
#include "tl_matrix_utils.h"

void showHelp(const std::string& progname) {
    std::cout << TlUtils::format("%s [options] MATRIX_FILE", progname.c_str()) << std::endl;
    std::cout << " OPTIONS:" << std::endl;
    std::cout << "  -h:      show help" << std::endl;
}

int main(int argc, char* argv[]) {
    TlGetopt opt(argc, argv, "h");

    if (opt["h"] == "defined") {
        showHelp(opt[0]);
        return EXIT_SUCCESS;
    }

    std::string path = opt[1];

    TlMatrixObject::HeaderInfo headerInfo;
    const int headerSize = TlMatrixUtils::getHeaderInfo(path, &headerInfo);
    if (headerSize > 0) {
        switch (headerInfo.matrixType) {
            case TlMatrixObject::RLHD:
                std::cout << "type: symmetric" << std::endl;
                std::cout << "row: " << headerInfo.numOfRows << std::endl;
                std::cout << "col: " << headerInfo.numOfCols << std::endl;
                break;

            case TlMatrixObject::CSFD:
                std::cout << "type: normal (column-major)" << std::endl;
                std::cout << "row: " << headerInfo.numOfRows << std::endl;
                std::cout << "col: " << headerInfo.numOfCols << std::endl;
                break;

            case TlMatrixObject::RSFD:
                std::cout << "type: normal (row-major)" << std::endl;
                std::cout << "row: " << headerInfo.numOfRows << std::endl;
                std::cout << "col: " << headerInfo.numOfCols << std::endl;
                break;

            case TlMatrixObject::ABGD:
                std::cout << "type: General Dens-matrix stored by Arrays Blocks" << std::endl;
                std::cout << "#vectors: " << headerInfo.numOfVectors << std::endl;
                std::cout << "sizeOfVector: " << headerInfo.sizeOfVector << std::endl;
                std::cout << TlUtils::format("unit: %d/%d", headerInfo.subunitId + 1, headerInfo.numOfSubunits)
                          << std::endl;
                std::cout << "chunk size: " << headerInfo.sizeOfChunk << std::endl;
                break;

            default:
                std::cerr << "unknown matrix type: " << path << std::endl;
        }
    } else {
        std::cerr << "can not open file: " << path << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
