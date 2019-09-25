#ifndef TOOL_COMMON_H
#define TOOL_COMMON_H

#include <stdio.h>
#include <time.h>
#include <cmath>
#include <cstdlib>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif  // _OPENMP

#ifdef HAVE_LAPACK
#include "tl_dense_general_matrix_lapack.h"
#include "tl_dense_symmetric_matrix_lapack.h"
#include "tl_dense_vector_lapack.h"
#endif  // HAVE_LAPACK

#ifdef HAVE_EIGEN
#include "tl_dense_general_matrix_eigen.h"
#include "tl_dense_symmetric_matrix_eigen.h"
#include "tl_dense_vector_eigen.h"
#endif  // HAVE_EIGEN

#ifdef HAVE_VIENNACL
#include "tl_dense_general_matrix_viennacl.h"
#include "tl_dense_symmetric_matrix_viennacl.h"
#include "tl_dense_vector_viennacl.h"

#include "tl_viennacl.h"
#include "viennacl/ocl/backend.hpp"
#include "viennacl/ocl/context.hpp"
#include "viennacl/ocl/enqueue.hpp"
#include "viennacl/ocl/platform.hpp"
#endif  // HAVE_VIENNACL

// -----------------------------------------------------------------------------
// General Matrix
// -----------------------------------------------------------------------------
template <typename DenseGeneralMatrix>
DenseGeneralMatrix getGeneralMatrixA(const TlMatrixObject::index_type row,
                                     const TlMatrixObject::index_type col) {
    const std::size_t numOfElements = row * col;
    DenseGeneralMatrix A(row, col);

    std::size_t i = 0;
    for (TlMatrixObject::index_type r = 0; r < row; ++r) {
        for (TlMatrixObject::index_type c = 0; c < col; ++c) {
            const double v = double(i) / double(numOfElements);
            A.set(r, c, v);
            ++i;
        }
    }

    return A;
};

template <typename DenseGeneralMatrix>
DenseGeneralMatrix getGeneralMatrixB(const TlMatrixObject::index_type row,
                                     const TlMatrixObject::index_type col) {
    const std::size_t numOfElements = row * col;
    DenseGeneralMatrix B(row, col);

    std::size_t i = numOfElements;
    for (TlMatrixObject::index_type r = 0; r < row; ++r) {
        for (TlMatrixObject::index_type c = 0; c < col; ++c) {
            const double v = double(i) / double(numOfElements);
            B.set(r, c, v);
            --i;
        }
    }

    return B;
};

#ifdef HAVE_VIENNACL
template <>
TlDenseGeneralMatrix_ViennaCL getGeneralMatrixA<TlDenseGeneralMatrix_ViennaCL>(
    const TlMatrixObject::index_type row,
    const TlMatrixObject::index_type col) {
    TlDenseGeneralMatrix_Eigen tmp =
        getGeneralMatrixA<TlDenseGeneralMatrix_Eigen>(row, col);
    TlDenseGeneralMatrix_ViennaCL answer(tmp);

    return answer;
};

template <>
TlDenseGeneralMatrix_ViennaCL getGeneralMatrixB<TlDenseGeneralMatrix_ViennaCL>(
    const TlMatrixObject::index_type row,
    const TlMatrixObject::index_type col) {
    TlDenseGeneralMatrix_Eigen tmp =
        getGeneralMatrixB<TlDenseGeneralMatrix_Eigen>(row, col);
    TlDenseGeneralMatrix_ViennaCL answer(tmp);

    return answer;
};
#endif  // HAVE_VIENNACL

// -----------------------------------------------------------------------------
template <typename DenseMatrixType>
void setupGeneralMatrix(const TlMatrixObject::index_type dim,
                        DenseMatrixType* pMat) {
    std::size_t count = 1;
    for (TlMatrixObject::index_type r = 0; r < dim; ++r) {
        for (TlMatrixObject::index_type c = 0; c < dim; ++c) {
            const double v = double(rand()) / double((RAND_MAX));
            pMat->set(r, c, v);
        }
    }
}

#ifdef HAVE_EIGEN
template <typename DenseMatrixType>
void setupGeneralMatrixWithEigen(const TlMatrixObject::index_type dim,
                                 DenseMatrixType* pMat) {
    TlDenseGeneralMatrix_Eigen tmp(dim, dim);
    for (TlMatrixObject::index_type r = 0; r < dim; ++r) {
        for (TlMatrixObject::index_type c = 0; c < dim; ++c) {
            const double v = double(rand()) / double((RAND_MAX));
            tmp.set(r, c, v);
        }
    }

    *pMat = tmp;
}
#endif  // HAVE_EIGEN

// -----------------------------------------------------------------------------
// Symmetric Matrix
// -----------------------------------------------------------------------------
template <typename SymmetricMatrixType>
void setupSymmetricMatrix(const TlMatrixObject::index_type dim,
                          SymmetricMatrixType* pMat) {
    for (TlMatrixObject::index_type r = 0; r < dim; ++r) {
        for (TlMatrixObject::index_type c = 0; c < r; ++c) {
            const double v = double(rand()) / double((RAND_MAX));
            pMat->set(r, c, v);
        }
        pMat->set(r, r, double(r));
    }
}

#ifdef HAVE_EIGEN
template <typename SymmetricMatrixType>
void setupSymmetricMatrixWithEigen(const TlMatrixObject::index_type dim,
                                   SymmetricMatrixType* pMat) {
    TlDenseSymmetricMatrix_Eigen tmp(dim);
    for (TlMatrixObject::index_type r = 0; r < dim; ++r) {
        for (TlMatrixObject::index_type c = 0; c < r; ++c) {
            const double v = double(rand()) / double((RAND_MAX));
            tmp.set(r, c, v);
        }
        tmp.set(r, r, double(r));
    }

    *pMat = tmp;
}
#endif  // HAVE_EIGEN

template <typename SymmetricMatrixType>
void setupSymmetricMatrixForDiagonal(const TlMatrixObject::index_type dim,
                                     SymmetricMatrixType* pMat) {
    pMat->resize(dim);
    for (TlMatrixObject::index_type r = 0; r < dim; ++r) {
        for (TlMatrixObject::index_type c = 0; c < r; ++c) {
            const double v = std::sqrt(0.3 * double(r) * double(c));
            pMat->set(r, c, v);
        }
        pMat->set(r, r, double(r));
    }
}

// -----------------------------------------------------------------------------
// prepare
// -----------------------------------------------------------------------------
void initRand() { srand((unsigned int)time(NULL)); }

#endif  // TOOL_COMMON_H
