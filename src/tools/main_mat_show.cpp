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

#include "TlGetopt.h"
#include "TlUtils.h"
#include "tl_dense_general_matrix_arrays_roworiented.h"
#include "tl_dense_general_matrix_lapack.h"
#include "tl_dense_symmetric_matrix_lapack.h"
#include "tl_matrix_utils.h"

int main(int argc, char* argv[]) {
    TlGetopt opt(argc, argv, "cghv");

    const bool bVerbose = (opt["v"] == "defined");
    const bool bCsvMode = (opt["c"] == "defined");
    const bool bGuessMode = (opt["g"] == "defined");

    std::string sPath = opt[1];
    if (bVerbose) {
        std::cerr << "loading... " << sPath << std::endl;
    }

    if (TlMatrixUtils::isLoadable(sPath, TlMatrixObject::RLHD) == true) {
        TlDenseSymmetricMatrix_Lapack M;
        M.load(sPath);

        if (bCsvMode == true) {
            M.saveCsv(std::cout);
        } else if (bGuessMode == true) {
            M.saveText(std::cout);
        } else {
            std::cout << M << std::endl;
        }
    } else if (TlMatrixUtils::isLoadable(sPath, TlMatrixObject::CSFD) == true) {
        TlDenseGeneralMatrix_Lapack M;
        M.load(sPath);

        if (bCsvMode == true) {
            M.saveCsv(std::cout);
        } else if (bGuessMode == true) {
            M.saveText(std::cout);
        } else {
            std::cout << M << std::endl;
        }
    } else if (TlMatrixUtils::isLoadable(sPath, TlMatrixObject::ABGD) == true) {
        TlMatrixObject::MatrixType type;
        TlMatrixObject::index_type numOfRows = 0;
        TlMatrixObject::index_type numOfCols = 0;
        std::size_t numOfItems = 0;
        TlMatrixObject::index_type numOfSubunits = 0;
        TlMatrixObject::index_type subunitId = 0;
        TlMatrixObject::index_type sizeOfChunk = 0;
        TlMatrixUtils::getHeaderInfo(sPath, &type, &numOfRows, &numOfCols,
                                     &numOfItems, &numOfSubunits, &subunitId,
                                     &sizeOfChunk);
        std::cout << "row: " << numOfRows << std::endl;
        std::cout << "col: " << numOfCols << std::endl;
        std::cout << TlUtils::format("unit: %d/%d", subunitId + 1,
                                     numOfSubunits)
                  << std::endl;
        std::cout << "chunk size: " << sizeOfChunk << std::endl;

        TlDenseGeneralMatrix_arrays_RowOriented M;
        bool isLoaded = M.load(sPath, -1);
        if (isLoaded) {
            std::cout << M << std::endl;
        } else {
            std::cout << "cannot load: " << sPath << std::endl;
        }

    } else {
        std::cerr << "unknown file type: " << sPath << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
