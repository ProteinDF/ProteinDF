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

#include <cassert>
#include <cstdlib>
#include <iostream>

#include "TlGetopt.h"
#include "TlUtils.h"
#include "tl_dense_general_matrix_lapack.h"
#include "tl_dense_symmetric_matrix_lapack.h"
#include "tl_matrix_object.h"
#include "tl_matrix_utils.h"

enum ExtendMode {
    EXTEND_UNDEFINED,
    EXTEND_ROW_WISE,
    EXTEND_COL_WISE,
    EXTEND_DIAGONAL_BLOCK
};

template <class OutputMatrixType>
int extend_diagonal_block_generalType(const TlMatrixObject* pMat1,
                                      const TlMatrixObject* pMat2,
                                      const std::string& outputMatrixPath,
                                      const bool isVerbose) {
    const TlMatrixObject::index_type row1 = pMat1->getNumOfRows();
    const TlMatrixObject::index_type col1 = pMat1->getNumOfCols();
    const TlMatrixObject::index_type row2 = pMat2->getNumOfRows();
    const TlMatrixObject::index_type col2 = pMat2->getNumOfCols();

    const TlMatrixObject::index_type newRow = row1 + row2;
    const TlMatrixObject::index_type newCol = col1 + col2;

    if ((newRow == 0) || (newCol == 0)) {
        std::cerr
            << TlUtils::format(
                   "need not resize matrix, because the size of new matrix "
                   "is (%d x %d)",
                   newRow, newCol)
            << std::endl;
        return EXIT_SUCCESS;
    }

    {
        OutputMatrixType mat3(newRow, newCol);

        // input1
        for (TlMatrixObject::index_type r = 0; r < row1; ++r) {
            for (TlMatrixObject::index_type c = 0; c < col1; ++c) {
                mat3.set(r, c, pMat1->get(r, c));
            }
        }

        // input2
        for (TlMatrixObject::index_type r = 0; r < row2; ++r) {
            for (TlMatrixObject::index_type c = 0; c < col2; ++c) {
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

template <class OutputMatrixType>
int extend_diagonal_block_symmetricType(const TlMatrixObject* pMat1,
                                        const TlMatrixObject* pMat2,
                                        const std::string& outputMatrixPath,
                                        const bool isVerbose) {
    const TlMatrixObject::index_type row1 = pMat1->getNumOfRows();
    const TlMatrixObject::index_type col1 = pMat1->getNumOfCols();
    const TlMatrixObject::index_type row2 = pMat2->getNumOfRows();
    const TlMatrixObject::index_type col2 = pMat2->getNumOfCols();

    const TlMatrixObject::index_type newRow = row1 + row2;
    const TlMatrixObject::index_type newCol = col1 + col2;

    if ((newRow == 0) || (newCol == 0)) {
        std::cerr
            << TlUtils::format(
                   "need not resize matrix, because the size of new matrix "
                   "is (%d x %d)",
                   newRow, newCol)
            << std::endl;
        return EXIT_SUCCESS;
    }

    {
        assert(row1 == col1);
        assert(row2 == col2);
        assert(newRow == newCol);

        OutputMatrixType mat3(newRow);

        // inpt1
        for (TlMatrixObject::index_type r = 0; r < row1; ++r) {
            for (TlMatrixObject::index_type c = 0; c <= r; ++c) {
                mat3.set(r, c, pMat1->get(r, c));
            }
        }

        // input2
        for (TlMatrixObject::index_type r = 0; r < row2; ++r) {
            for (TlMatrixObject::index_type c = 0; c <= r; ++c) {
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

template <class OutputMatrixType>
int extend_row_wise(const TlMatrixObject* pMat1, const TlMatrixObject* pMat2,
                    const std::string& outputMatrixPath, const bool isVerbose) {
    const TlMatrixObject::index_type row1 = pMat1->getNumOfRows();
    const TlMatrixObject::index_type col1 = pMat1->getNumOfCols();
    const TlMatrixObject::index_type row2 = pMat2->getNumOfRows();
    const TlMatrixObject::index_type col2 = pMat2->getNumOfCols();

    if (col1 != col2) {
        std::cerr << TlUtils::format("The two col dims are mismatch. %d != %d",
                                     col1, col2)
                  << std::endl;
        return EXIT_FAILURE;
    }

    const TlMatrixObject::index_type newRow = row1 + row2;
    const TlMatrixObject::index_type newCol = col1;
    assert(col1 == col2);

    if ((newRow == 0) || (newCol == 0)) {
        std::cerr
            << TlUtils::format(
                   "need not resize matrix, because the size of new matrix "
                   "is (%d x %d)",
                   newRow, newCol)
            << std::endl;
        return EXIT_SUCCESS;
    }

    OutputMatrixType mat3(newRow, newCol);
    // input1
    for (TlMatrixObject::index_type r = 0; r < row1; ++r) {
        for (TlMatrixObject::index_type c = 0; c < col1; ++c) {
            const double v = pMat1->get(r, c);
            mat3.set(r, c, v);
        }
    }

    // input2
    for (TlMatrixObject::index_type r = 0; r < row2; ++r) {
        for (TlMatrixObject::index_type c = 0; c < newCol; ++c) {
            const double v = pMat2->get(r, c);
            mat3.set(row1 + r, c, v);
        }
    }

    mat3.save(outputMatrixPath);
    return EXIT_SUCCESS;
}

template <class OutputMatrixType>
int extend_column_wise(const TlMatrixObject* pMat1, const TlMatrixObject* pMat2,
                       const std::string& outputMatrixPath,
                       const bool isVerbose) {
    const TlMatrixObject::index_type row1 = pMat1->getNumOfRows();
    const TlMatrixObject::index_type col1 = pMat1->getNumOfCols();
    const TlMatrixObject::index_type row2 = pMat2->getNumOfRows();
    const TlMatrixObject::index_type col2 = pMat2->getNumOfCols();

    if (row1 != row2) {
        std::cerr << TlUtils::format("The two row dims are mismatch. %d != %d",
                                     row1, row2)
                  << std::endl;
        return EXIT_FAILURE;
    }

    const TlMatrixObject::index_type newRow = row1;
    assert(row1 == row2);
    const TlMatrixObject::index_type newCol = col1 + col2;

    if ((newRow == 0) || (newCol == 0)) {
        std::cerr
            << TlUtils::format(
                   "need not resize matrix, because the size of new matrix "
                   "is (%d x %d)",
                   newRow, newCol)
            << std::endl;
        return EXIT_SUCCESS;
    }

    OutputMatrixType mat3(newRow, newCol);

    // input1
    for (TlMatrixObject::index_type r = 0; r < row1; ++r) {
        for (TlMatrixObject::index_type c = 0; c < col1; ++c) {
            const double v = pMat1->get(r, c);
            mat3.set(r, c, v);
        }
    }

    // input2
    for (TlMatrixObject::index_type c = 0; c < col2; ++c) {
        for (TlMatrixObject::index_type r = 0; r < newRow; ++r) {
            const double v = pMat2->get(r, c);
            mat3.set(r, col1 + c, v);
        }
    }

    mat3.save(outputMatrixPath);
    return EXIT_SUCCESS;
}

// -----------------------------------------------------------------------------
// Main
// -----------------------------------------------------------------------------
void showHelp(const std::string& name) {
    std::cout << TlUtils::format(
                     "%s extend mode [options] base_matrix_path "
                     "reference_matrix_path  output_path",
                     name.c_str())
              << std::endl;
    std::cout << " OPTIONS:" << std::endl;
    std::cout << "extend mode" << std::endl;
    std::cout << "  -d:      extend diagonal block" << std::endl;
    std::cout << "  -r:      extend row-wise" << std::endl;
    std::cout << "  -c:      extend column-wise" << std::endl;
    std::cout << std::endl;
    std::cout << "  -h:      show help" << std::endl;
    std::cout << "  -v:      verbose" << std::endl;
}

int main(int argc, char* argv[]) {
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

    TlMatrixObject* pMat1 = NULL;
    bool isSymmetric1 = false;
    if (isVerbose == true) {
        std::cerr << "load matrix: " << baseMatrixPath << std::endl;
    }
    if (TlMatrixUtils::isLoadable(baseMatrixPath, TlMatrixObject::RLHD)) {
        // TlDenseSymmetricMatrix_Lapack tmp;
        // tmp.load(baseMatrixPath);
        // pMat1 = new TlDenseGeneralMatrix_Lapack(tmp);
        pMat1 = new TlDenseSymmetricMatrix_Lapack();
        pMat1->load(baseMatrixPath);
        isSymmetric1 = true;
    } else if (TlMatrixUtils::isLoadable(baseMatrixPath,
                                         TlMatrixObject::CSFD)) {
        pMat1 = new TlDenseGeneralMatrix_Lapack();
        pMat1->load(baseMatrixPath);
    } else {
        std::cerr << TlUtils::format("cannot load: %s", baseMatrixPath.c_str())
                  << std::endl;
        std::cerr << "create new matrix." << std::endl;

        copyMode = true;
    }

    TlMatrixObject* pMat2 = NULL;
    bool isSymmetric2 = false;
    if (isVerbose == true) {
        std::cerr << "load matrix: " << refMatrixPath << std::endl;
    }
    if (TlMatrixUtils::isLoadable(refMatrixPath, TlMatrixObject::RLHD)) {
        // TlDenseSymmetricMatrix_Lapack tmp;
        // tmp.load(refMatrixPath);
        pMat2 = new TlDenseSymmetricMatrix_Lapack();
        pMat2->load(refMatrixPath);
        isSymmetric2 = true;
    } else if (TlMatrixUtils::isLoadable(refMatrixPath, TlMatrixObject::CSFD)) {
        pMat2 = new TlDenseGeneralMatrix_Lapack();
        pMat2->load(refMatrixPath);
    } else {
        std::cerr << TlUtils::format("cannot load: %s", refMatrixPath.c_str())
                  << std::endl;
        return EXIT_FAILURE;
    }

    int retval = EXIT_FAILURE;
    if (copyMode) {
        pMat2->save(outputMatrixPath);
        retval = EXIT_SUCCESS;
    } else {
        switch (extendMode) {
            case EXTEND_ROW_WISE:
                retval = extend_row_wise<TlDenseGeneralMatrix_Lapack>(pMat1, pMat2, outputMatrixPath, isVerbose);
                break;

            case EXTEND_COL_WISE:
                retval = extend_column_wise<TlDenseGeneralMatrix_Lapack>(pMat1, pMat2, outputMatrixPath, isVerbose);
                break;

            case EXTEND_DIAGONAL_BLOCK: {
                if ((isSymmetric1 == true) && (isSymmetric2 == true)) {
                    retval = extend_diagonal_block_symmetricType<TlDenseSymmetricMatrix_Lapack>(pMat1, pMat2, outputMatrixPath, isVerbose);
                } else {
                    retval = extend_diagonal_block_generalType<TlDenseGeneralMatrix_Lapack>(pMat1, pMat2, outputMatrixPath, isVerbose);
                }
            } break;

            default:
                std::cerr << "program error. " << __FILE__ << __LINE__
                          << std::endl;
                break;
        }
    }

    delete pMat1;
    pMat1 = NULL;
    delete pMat2;
    pMat2 = NULL;

    return retval;
}
