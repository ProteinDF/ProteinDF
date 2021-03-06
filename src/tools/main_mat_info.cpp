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
    std::cout << TlUtils::format("%s [options] MATRIX_FILE", progname.c_str())
              << std::endl;
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

    TlMatrixObject::MatrixType type;
    TlMatrixObject::index_type numOfRows = 0;
    TlMatrixObject::index_type numOfCols = 0;

    if (TlMatrixUtils::isLoadable(path, TlMatrixObject::RLHD) == true) {
        TlMatrixUtils::getHeaderInfo(path, &type, &numOfRows, &numOfCols);

        std::cout << "type: symmetric" << std::endl;
        std::cout << "row: " << numOfRows << std::endl;
        std::cout << "col: " << numOfCols << std::endl;
    } else if (TlMatrixUtils::isLoadable(path, TlMatrixObject::CSFD) == true) {
        TlMatrixUtils::getHeaderInfo(path, &type, &numOfRows, &numOfCols);

        std::cout << "type: normal (column-major)" << std::endl;
        std::cout << "row: " << numOfRows << std::endl;
        std::cout << "col: " << numOfCols << std::endl;
    } else if (TlMatrixUtils::isLoadable(path, TlMatrixObject::RSFD) == true) {
        TlMatrixUtils::getHeaderInfo(path, &type, &numOfRows, &numOfCols);

        std::cout << "type: normal (row-major)" << std::endl;
        std::cout << "row: " << numOfRows << std::endl;
        std::cout << "col: " << numOfCols << std::endl;
    } else if (TlMatrixUtils::isLoadable(path, TlMatrixObject::ABGD) == true) {
        std::size_t numOfItems = 0;
        TlMatrixObject::index_type numOfSubunits = 0;
        TlMatrixObject::index_type subunitId = 0;
        TlMatrixObject::index_type sizeOfChunk = 0;
        TlMatrixUtils::getHeaderInfo(path, &type, &numOfRows, &numOfCols,
                                     &numOfItems, &numOfSubunits, &subunitId,
                                     &sizeOfChunk);
        std::cout << "type: General Dens-matrix stored by Arrays Blocks"
                  << std::endl;
        std::cout << "row: " << numOfRows << std::endl;
        std::cout << "col: " << numOfCols << std::endl;
        std::cout << TlUtils::format("unit: %d/%d", subunitId + 1,
                                     numOfSubunits)
                  << std::endl;
        std::cout << "chunk size: " << sizeOfChunk << std::endl;
    } else {
        std::cerr << "can not open file: " << path << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
