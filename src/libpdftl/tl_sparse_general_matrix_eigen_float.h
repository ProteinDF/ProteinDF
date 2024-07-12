#ifndef TL_SPARSE_GENERAL_MATRIX_EIGEN_FLOAT_H
#define TL_SPARSE_GENERAL_MATRIX_EIGEN_FLOAT_H

#include "tl_dense_general_matrix_eigen_float.h"
#include "tl_dense_vector_eigen_float.h"
#include "tl_sparse_general_matrix_impl_eigen_float.h"
#include "tl_sparse_general_matrix_object.h"

class TlSparseSymmetricMatrix_EigenFloat;
class TlSparseSymmetricMatrix_ImplEigenFloat;
class TlSparseGeneralMatrix_ViennaCL;

class TlSparseGeneralMatrix_EigenFloat : public TlSparseGeneralMatrixObject {
    // ---------------------------------------------------------------------------
    // constructor & destructor
    // ---------------------------------------------------------------------------
public:
    explicit TlSparseGeneralMatrix_EigenFloat(
        const TlMatrixObject::index_type row = 1,
        const TlMatrixObject::index_type col = 1);
    TlSparseGeneralMatrix_EigenFloat(const TlSparseGeneralMatrix_EigenFloat& rhs);
    TlSparseGeneralMatrix_EigenFloat(const TlSparseSymmetricMatrix_EigenFloat& rhs);
    TlSparseGeneralMatrix_EigenFloat(const TlSparseGeneralMatrix_ImplEigenFloat& rhs);
    TlSparseGeneralMatrix_EigenFloat(const TlDenseGeneralMatrix_EigenFloat& rhs);

#ifdef HAVE_VIENNACL
    TlSparseGeneralMatrix_EigenFloat(const TlSparseGeneralMatrix_ViennaCL& rhs);
#endif  // HAVE_VIENNACL

    virtual ~TlSparseGeneralMatrix_EigenFloat();

    // ---------------------------------------------------------------------------
    // operators
    // ---------------------------------------------------------------------------
public:
    TlSparseGeneralMatrix_EigenFloat& operator=(
        const TlSparseGeneralMatrix_EigenFloat& rhs);
    TlSparseGeneralMatrix_EigenFloat& operator+=(
        const TlSparseGeneralMatrix_EigenFloat& rhs);
    TlSparseGeneralMatrix_EigenFloat& operator-=(
        const TlSparseGeneralMatrix_EigenFloat& rhs);
    TlSparseGeneralMatrix_EigenFloat& operator*=(const double coef);

    // ---------------------------------------------------------------------------
    // others
    // ---------------------------------------------------------------------------
    friend class TlDenseGeneralMatrix_EigenFloat;
    friend class TlSparseSymmetricMatrix_EigenFloat;
    friend class TlSparseSymmetricMatrix_ImplEigenFloat;
    friend class TlSparseGeneralMatrix_ViennaCL;

    // DM(G) = DM(G) * SM(G)
    friend TlDenseGeneralMatrix_EigenFloat operator*(
        const TlDenseGeneralMatrix_EigenFloat& mat1,
        const TlSparseGeneralMatrix_EigenFloat& mat2);
    // DM(G) = SM(G) * DM(G)
    friend TlDenseGeneralMatrix_EigenFloat operator*(
        const TlSparseGeneralMatrix_EigenFloat& mat1,
        const TlDenseGeneralMatrix_EigenFloat& mat2);

    // SM(G) = SM(G) * SM(S)
    friend TlSparseGeneralMatrix_EigenFloat operator*(
        const TlSparseGeneralMatrix_EigenFloat& sm1,
        const TlSparseSymmetricMatrix_EigenFloat& sm2);
    // SM(G) = SM(S) * SM(G)
    friend TlSparseGeneralMatrix_EigenFloat operator*(
        const TlSparseSymmetricMatrix_EigenFloat& sm1,
        const TlSparseGeneralMatrix_EigenFloat& sm2);
    // SM(G) = SM(G) * SM(G)
    friend TlSparseGeneralMatrix_EigenFloat operator*(
        const TlSparseGeneralMatrix_EigenFloat& sm1,
        const TlSparseGeneralMatrix_EigenFloat& sm2);

    // DV = SM(G) * DV
    friend TlDenseVector_EigenFloat operator*(const TlSparseGeneralMatrix_EigenFloat& mat,
                                              const TlDenseVector_EigenFloat& vtr);
};
#endif  // TL_SPARSE_GENERAL_MATRIX_EIGEN_FLOAT_H
