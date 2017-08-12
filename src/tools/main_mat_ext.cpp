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

enum ExtendMode {
    EXTEND_UNDEFINED,
    EXTEND_ROW_WISE,
    EXTEND_COL_WISE,
    EXTEND_DIAGONAL_BLOCK
};


int extend_diagonal_block(const TlMatrix* pMat1, const TlMatrix* pMat2,
                          const std::string& outputMatrixPath,
                          const bool isSaveSymmetric,
                          const bool isVerbose);
int extend_row_wise(const TlMatrix* pMat1, const TlMatrix* pMat2,
                    const std::string& outputMatrixPath,
                    const bool isVerbose);
int extend_column_wise(const TlMatrix* pMat1, const TlMatrix* pMat2,
                       const std::string& outputMatrixPath,
                       const bool isVerbose);
void showHelp(const std::string& name);

int main(int argc, char* argv[])
{
    TlGetopt opt(argc, argv, "drchv");
    
    if (opt["h"] == "defined") {
        showHelp(opt[0]);
        return EXIT_SUCCESS;
    }
    
    const bool isVerbose = (opt["v"] == "defined");
    if (opt.getCount() <= 2) {
        showHelp(opt[0]);
        return EXIT_FAILURE;
    }
    std::string baseMatrixPath = opt[1];
    std::string refMatrixPath = opt[2];
    std::string outputMatrixPath = opt[3];

    ExtendMode extendMode = EXTEND_UNDEFINED;
    if (opt["d"] == "defined") {
        extendMode = EXTEND_DIAGONAL_BLOCK;
    } else if (opt["r"] == "defined") {
        extendMode = EXTEND_ROW_WISE;
    } else if (opt["c"] == "defined") {
        extendMode = EXTEND_COL_WISE;
    }
    if (extendMode == EXTEND_UNDEFINED) {
        std::cerr << "please set extend mode (-d, -r, or -c)." << std::endl;
        showHelp(opt[0]);
        return EXIT_FAILURE;
    }


    bool copyMode = false;

    TlMatrix* pMat1 = NULL;
    TlMatrix::index_type row1 = 0;
    TlMatrix::index_type col1 = 0;
    bool isSymmetric1 = false;
    if (isVerbose == true) {
        std::cerr << "load matrix: " << baseMatrixPath << std::endl;
    }
    if (TlSymmetricMatrix::isLoadable(baseMatrixPath)) {
        pMat1 = new TlSymmetricMatrix();
        isSymmetric1 = true;
        
        pMat1->load(baseMatrixPath);
        row1 = pMat1->getNumOfRows();
        col1 = pMat1->getNumOfCols();
    } else if (TlMatrix::isLoadable(baseMatrixPath)) {
        pMat1 = new TlMatrix();
        
        pMat1->load(baseMatrixPath);
        row1 = pMat1->getNumOfRows();
        col1 = pMat1->getNumOfCols();
    } else {
        std::cerr << TlUtils::format("cannot load: %s", baseMatrixPath.c_str()) << std::endl;
        std::cerr << "create new matrix." << std::endl;
        
        copyMode = true;
    }
    
    
    TlMatrix* pMat2 = NULL;
    TlMatrix::index_type row2 = 0;
    TlMatrix::index_type col2 = 0;
    bool isSymmetric2 = false;
    if (isVerbose == true) {
        std::cerr << "load matrix: " << refMatrixPath << std::endl;
    }
    if (TlSymmetricMatrix::isLoadable(refMatrixPath)) {
        pMat2 = new TlSymmetricMatrix();
        isSymmetric2 = true;
    } else if (TlMatrix::isLoadable(refMatrixPath)) {
        pMat2 = new TlMatrix();
    } else {
        std::cerr << TlUtils::format("cannot load: %s", refMatrixPath.c_str()) << std::endl;
        return EXIT_FAILURE;
    }
    pMat2->load(refMatrixPath);
    row2 = pMat2->getNumOfRows();
    col2 = pMat2->getNumOfCols();

    int retval = EXIT_FAILURE;
    if (copyMode) {
        pMat2->save(outputMatrixPath);
        retval = EXIT_SUCCESS;
    } else {
        switch (extendMode)
        {
        case EXTEND_ROW_WISE:
            retval = extend_row_wise(pMat1, pMat2, outputMatrixPath, isVerbose);
            break;

        case EXTEND_COL_WISE:
            retval = extend_column_wise(pMat1, pMat2, outputMatrixPath, isVerbose);
            break;

        case EXTEND_DIAGONAL_BLOCK:
            retval = extend_diagonal_block(pMat1, pMat2, outputMatrixPath,
                                           ((isSymmetric1 == true) && (isSymmetric2 == true)),
                                           isVerbose);
            break;

        default:
            std::cerr << "program error. " << __FILE__ << __LINE__ << std::endl;
            break;
        }
    }
    
    delete pMat1;
    pMat1 = NULL;
    delete pMat2;
    pMat2 = NULL;
    
    return EXIT_SUCCESS;
}


int extend_diagonal_block(const TlMatrix* pMat1, const TlMatrix* pMat2,
                          const std::string& outputMatrixPath,
                          const bool isSaveSymmetric,
                          const bool isVerbose)
{
    const TlMatrix::index_type row1 = pMat1->getNumOfRows();
    const TlMatrix::index_type col1 = pMat1->getNumOfCols();
    const TlMatrix::index_type row2 = pMat2->getNumOfRows();
    const TlMatrix::index_type col2 = pMat2->getNumOfCols();

    const TlMatrix::index_type newRow = row1 + row2;
    const TlMatrix::index_type newCol = col1 + col2;

    if ((newRow == 0) || (newCol == 0)) {
        std::cerr << TlUtils::format("need not resize matrix, because the size of new matrix is (%d x %d)", newRow, newCol)
                  << std::endl;
        return EXIT_SUCCESS;
    }
    
    if (isSaveSymmetric) {
        assert(row1 == col1);
        assert(row2 == col2);
        assert(newRow == newCol);
        
        TlSymmetricMatrix mat3 = *pMat1;
        mat3.resize(newRow);
        
        for (TlMatrix::index_type r = 0; r < row2; ++r) {
            for (TlMatrix::index_type c = 0; c <= r; ++c) {
                mat3.set(row1 + r, col1 + c, pMat2->get(r, c));
            }
        }
        
        if (isVerbose == true) {
            std::cerr << "save matrix: " << outputMatrixPath << std::endl;
        }
        mat3.save(outputMatrixPath);
    } else {
        TlMatrix mat3 = *pMat1;
        mat3.resize(newRow, newCol);
        
        for (TlMatrix::index_type r = 0; r < row2; ++r) {
            for (TlMatrix::index_type c = 0; c < col2; ++c) {
                mat3.set(row1 + r, col1 + c, pMat2->get(r, c));
            }
        }
        
        if (isVerbose == true) {
            std::cerr << "save matrix: " << outputMatrixPath << std::endl;
        }
        mat3.save(outputMatrixPath);
    }

    return EXIT_SUCCESS;
}


int extend_row_wise(const TlMatrix* pMat1, const TlMatrix* pMat2,
                    const std::string& outputMatrixPath,
                    const bool isVerbose)
{
    const TlMatrix::index_type row1 = pMat1->getNumOfRows();
    const TlMatrix::index_type col1 = pMat1->getNumOfCols();
    const TlMatrix::index_type row2 = pMat2->getNumOfRows();
    const TlMatrix::index_type col2 = pMat2->getNumOfCols();

    if (col1 != col2) {
        std::cerr << TlUtils::format("The two col dims are mismatch. %d != %d",
                                     col1, col2) 
                  << std::endl;
        return EXIT_FAILURE;
    }
    
    const TlMatrix::index_type newRow = row1 + row2;
    const TlMatrix::index_type newCol = col1;
    assert(col1 == col2);
    
    if ((newRow == 0) || (newCol == 0)) {
        std::cerr << TlUtils::format("need not resize matrix, because the size of new matrix is (%d x %d)", newRow, newCol)
                  << std::endl;
        return EXIT_SUCCESS;
    }
    
    TlMatrix mat3 = *pMat1;
    mat3.resize(newRow, newCol);
    for (TlMatrix::index_type r = 0; r < row2; ++r) {
        for (TlMatrix::index_type c = 0; c < newCol; ++c) {
            const double v = pMat2->get(r, c);
            mat3.set(row1 + r, c, v);
        }
    }
    
    mat3.save(outputMatrixPath);
    return EXIT_SUCCESS;
}

int extend_column_wise(const TlMatrix* pMat1, const TlMatrix* pMat2,
                       const std::string& outputMatrixPath,
                       const bool isVerbose)
{
    const TlMatrix::index_type row1 = pMat1->getNumOfRows();
    const TlMatrix::index_type col1 = pMat1->getNumOfCols();
    const TlMatrix::index_type row2 = pMat2->getNumOfRows();
    const TlMatrix::index_type col2 = pMat2->getNumOfCols();

    if (row1 != row2) {
        std::cerr << TlUtils::format("The two row dims are mismatch. %d != %d",
                                     row1, row2) 
                  << std::endl;
        return EXIT_FAILURE;
    }
    
    const TlMatrix::index_type newRow = row1;
    assert(row1 == row2);
    const TlMatrix::index_type newCol = col1 + col2;
    
    if ((newRow == 0) || (newCol == 0)) {
        std::cerr << TlUtils::format("need not resize matrix, because the size of new matrix is (%d x %d)", newRow, newCol)
                  << std::endl;
        return EXIT_SUCCESS;
    }
    
    TlMatrix mat3 = *pMat1;
    mat3.resize(newRow, newCol);
    for (TlMatrix::index_type c = 0; c < col2; ++c) {
        for (TlMatrix::index_type r = 0; r < newRow; ++r) {
            const double v = pMat2->get(r, c);
            mat3.set(r, col1 + c, v);
        }
    }
    
    mat3.save(outputMatrixPath);
    return EXIT_SUCCESS;
}

void showHelp(const std::string& name)
{
    std::cout << TlUtils::format("%s extend mode [options] base_matrix_path reference_matrix_path  output_path", name.c_str()) << std::endl;
    std::cout << " OPTIONS:" << std::endl;
    std::cout << "extend mode" << std::endl;
    std::cout << "  -d:      extend diagonal block" << std::endl;
    std::cout << "  -r:      extend row-wise" << std::endl;
    std::cout << "  -c:      extend coloumn-wise" << std::endl;
    std::cout << std::endl;
    std::cout << "  -h:      show help" << std::endl;
    std::cout << "  -v:      verbose" << std::endl;
}

