#ifndef TL_SPARSE_GENERAL_MATRIX_VIENNACL_FLOAT_H
#define TL_SPARSE_GENERAL_MATRIX_VIENNACL_FLOAT_H

#include "tl_sparse_general_matrix_impl_viennacl.h"
#include "tl_sparse_general_matrix_object.h"

class TlSparseSymmetricMatrix_ViennaCLFloat;
class TlSparseSymmetricMatrix_ImplViennaCLFloat;
class TlDenseGeneralMatrix_ViennaCLFloat;
class TlDenseVector_ViennaCLFloat;
class TlSparseGeneralMatrix_EigenFloat;

class TlSparseGeneralMatrix_ViennaCLFloat : public TlSparseGeneralMatrixObject {
    // ---------------------------------------------------------------------------
    // constructor & destructor
    // ---------------------------------------------------------------------------
public:
    explicit TlSparseGeneralMatrix_ViennaCLFloat(
        const TlMatrixObject::index_type row = 1,
        const TlMatrixObject::index_type col = 1);
    TlSparseGeneralMatrix_ViennaCLFloat(const TlSparseGeneralMatrix_ViennaCLFloat& rhs);
    TlSparseGeneralMatrix_ViennaCLFloat(const TlSparseSymmetricMatrix_ViennaCLFloat& rhs);
    TlSparseGeneralMatrix_ViennaCLFloat(const TlSparseGeneralMatrix_ImplViennaCLFloat& rhs);
    TlSparseGeneralMatrix_ViennaCLFloat(const TlDenseGeneralMatrix_ViennaCLFloat& rhs);

#ifdef HAVE_EIGEN
    TlSparseGeneralMatrix_ViennaCLFloat(const TlSparseGeneralMatrix_EigenFloat& rhs);
#endif  // HAVE_EIGEN

    virtual ~TlSparseGeneralMatrix_ViennaCLFloat();

    // ---------------------------------------------------------------------------
    // operators
    // ---------------------------------------------------------------------------
public:
    TlSparseGeneralMatrix_ViennaCLFloat& operator=(
        const TlSparseGeneralMatrix_ViennaCLFloat& rhs);
    TlSparseGeneralMatrix_ViennaCLFloat& operator+=(
        const TlSparseGeneralMatrix_ViennaCLFloat& rhs);
    TlSparseGeneralMatrix_ViennaCLFloat& operator-=(
        const TlSparseGeneralMatrix_ViennaCLFloat& rhs);
    TlSparseGeneralMatrix_ViennaCLFloat& operator*=(const double coef);

    // ---------------------------------------------------------------------------
    // others
    // ---------------------------------------------------------------------------
    friend class TlSparseSymmetricMatrix_ViennaCLFloat;
    friend class TlSparseSymmetricMatrix_ImplViennaCLFloat;
    friend class TlSparseGeneralMatrix_EigenFloat;
    friend class TlDenseGeneralMatrix_ViennaCLFloat;

    // SM(G) = SM(G) * SM(G)
    friend TlSparseGeneralMatrix_ViennaCLFloat operator*(
        const TlSparseGeneralMatrix_ViennaCLFloat& sm1,
        const TlSparseGeneralMatrix_ViennaCLFloat& sm2);
    // SM(G) = SM(G) * SM(S)
    friend TlSparseGeneralMatrix_ViennaCLFloat operator*(
        const TlSparseGeneralMatrix_ViennaCLFloat& sm1,
        const TlSparseSymmetricMatrix_ViennaCLFloat& sm2);
    // SM(G) = SM(S) * SM(G)
    friend TlSparseGeneralMatrix_ViennaCLFloat operator*(
        const TlSparseSymmetricMatrix_ViennaCLFloat& sm1,
        const TlSparseGeneralMatrix_ViennaCLFloat& sm2);
    // DV = SM(G) * DV
    friend TlDenseVector_ViennaCLFloat operator*(
        const TlSparseGeneralMatrix_ViennaCLFloat& mat,
        const TlDenseVector_ViennaCLFloat& vtr);
};
#endif  // TL_SPARSE_GENERAL_MATRIX_VIENNACL_FLOAT_H
