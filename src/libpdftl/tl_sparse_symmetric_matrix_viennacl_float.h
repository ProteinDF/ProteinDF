#ifndef TL_SPARSE_SYMMETRIC_MATRIX_VIENNACL_FLOAT_H
#define TL_SPARSE_SYMMETRIC_MATRIX_VIENNACL_FLOAT_H

#include "tl_sparse_symmetric_matrix_object.h"

class TlSparseGeneralMatrix_ViennaCLFloat;
class TlSparseGeneralMatrix_ImplViennaCLFloat;
class TlSparseSymmetricMatrix_EigenFloat;
class TlDenseSymmetricMatrix_ViennaCLFloat;

class TlSparseSymmetricMatrix_ViennaCLFloat : public TlSparseSymmetricMatrixObject {
    // ---------------------------------------------------------------------------
    // constructor & destructor
    // ---------------------------------------------------------------------------
public:
    explicit TlSparseSymmetricMatrix_ViennaCLFloat(const TlMatrixObject::index_type dim = 1);
    TlSparseSymmetricMatrix_ViennaCLFloat(const TlSparseSymmetricMatrix_ViennaCLFloat& rhs);
    TlSparseSymmetricMatrix_ViennaCLFloat(const TlSparseGeneralMatrix_ViennaCLFloat& rhs);
    TlSparseSymmetricMatrix_ViennaCLFloat(const TlDenseSymmetricMatrix_ViennaCLFloat& rhs);

#ifdef HAVE_EIGEN
    TlSparseSymmetricMatrix_ViennaCLFloat(const TlSparseSymmetricMatrix_EigenFloat& rhs);
#endif  // HAVE_EIGEN

    virtual ~TlSparseSymmetricMatrix_ViennaCLFloat();

    // ---------------------------------------------------------------------------
    // operators
    // ---------------------------------------------------------------------------
public:
    TlSparseSymmetricMatrix_ViennaCLFloat& operator=(
        const TlSparseSymmetricMatrix_ViennaCLFloat& rhs);
    TlSparseSymmetricMatrix_ViennaCLFloat& operator+=(
        const TlSparseSymmetricMatrix_ViennaCLFloat& rhs);
    TlSparseSymmetricMatrix_ViennaCLFloat& operator-=(
        const TlSparseSymmetricMatrix_ViennaCLFloat& rhs);
    TlSparseSymmetricMatrix_ViennaCLFloat& operator*=(const double coef);

    // ---------------------------------------------------------------------------
    // others
    // ---------------------------------------------------------------------------
    friend class TlSparseGeneralMatrix_ViennaCLFloat;
    friend class TlSparseGeneralMatrix_ImplViennaCLFloat;
    friend class TlSparseSymmetricMatrix_EigenFloat;
    friend class TlDenseSymmetricMatrix_ViennaCLFloat;

    friend TlSparseGeneralMatrix_ViennaCLFloat operator*(
        const TlSparseGeneralMatrix_ViennaCLFloat& sm1,
        const TlSparseSymmetricMatrix_ViennaCLFloat& sm2);
    friend TlSparseGeneralMatrix_ViennaCLFloat operator*(
        const TlSparseSymmetricMatrix_ViennaCLFloat& sm1,
        const TlSparseGeneralMatrix_ViennaCLFloat& sm2);
    friend TlSparseGeneralMatrix_ViennaCLFloat operator*(
        const TlSparseSymmetricMatrix_ViennaCLFloat& sm1,
        const TlSparseSymmetricMatrix_ViennaCLFloat& sm2);
};

#endif  // TL_SPARSE_SYMMETRIC_MATRIX_VIENNACL_FLOAT_H
