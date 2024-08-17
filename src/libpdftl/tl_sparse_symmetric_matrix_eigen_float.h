#ifndef TL_SPARSE_SYMMETRIC_MATRIX_EIGEN_FLOAT_H
#define TL_SPARSE_SYMMETRIC_MATRIX_EIGEN_FLOAT_H

#include "tl_sparse_symmetric_matrix_object.h"

class TlDenseGeneralMatrix_EigenFloat;
class TlDenseSymmetricMatrix_EigenFloat;
class TlDenseVector_EigenFloat;
class TlSparseGeneralMatrix_EigenFloat;
class TlSparseGeneralMatrix_ImplEigenFloat;
class TlSparseSymmetricMatrix_ViennaCLFloat;

class TlSparseSymmetricMatrix_EigenFloat : public TlSparseSymmetricMatrixObject {
    // ---------------------------------------------------------------------------
    // constructor & destructor
    // ---------------------------------------------------------------------------
public:
    explicit TlSparseSymmetricMatrix_EigenFloat(const TlMatrixObject::index_type dim = 1);
    TlSparseSymmetricMatrix_EigenFloat(const TlSparseSymmetricMatrix_EigenFloat& rhs);
    TlSparseSymmetricMatrix_EigenFloat(const TlSparseGeneralMatrix_EigenFloat& rhs);
    TlSparseSymmetricMatrix_EigenFloat(const TlDenseSymmetricMatrix_EigenFloat& rhs);

#ifdef HAVE_VIENNACL
    TlSparseSymmetricMatrix_EigenFloat(const TlSparseSymmetricMatrix_ViennaCLFloat& rhs);
#endif  // HAVE_VIENNACL

    virtual ~TlSparseSymmetricMatrix_EigenFloat();

    // ---------------------------------------------------------------------------
    // operators
    // ---------------------------------------------------------------------------
public:
    TlSparseSymmetricMatrix_EigenFloat& operator=(const TlSparseSymmetricMatrix_EigenFloat& rhs);

    TlSparseSymmetricMatrix_EigenFloat& operator+=(const TlSparseSymmetricMatrix_EigenFloat& sm);
    TlSparseSymmetricMatrix_EigenFloat& operator-=(const TlSparseSymmetricMatrix_EigenFloat& sm);

    // ---------------------------------------------------------------------------
    // others
    // ---------------------------------------------------------------------------
    friend class TlDenseSymmetricMatrix_EigenFloat;
    friend class TlSparseGeneralMatrix_EigenFloat;
    friend class TlSparseGeneralMatrix_ImplEigenFloat;
    friend class TlSparseSymmetricMatrix_ViennaCLFloat;

    // DM(G) = DM(G) * SM(S)
    friend TlDenseGeneralMatrix_EigenFloat operator*(const TlDenseGeneralMatrix_EigenFloat& mat1,
                                                     const TlSparseSymmetricMatrix_EigenFloat& mat2);
    // DM(G) = SM(S) * DM(G)
    friend TlDenseGeneralMatrix_EigenFloat operator*(const TlSparseSymmetricMatrix_EigenFloat& mat1,
                                                     const TlDenseGeneralMatrix_EigenFloat& mat2);

    // SM(G) = SM(G) * SM(S)
    friend TlSparseGeneralMatrix_EigenFloat operator*(const TlSparseGeneralMatrix_EigenFloat& sm1,
                                                      const TlSparseSymmetricMatrix_EigenFloat& sm2);
    // SM(G) = SM(S) * SM(G)
    friend TlSparseGeneralMatrix_EigenFloat operator*(const TlSparseSymmetricMatrix_EigenFloat& sm1,
                                                      const TlSparseGeneralMatrix_EigenFloat& sm2);
    // SM(G) = SM(S) * SM(S)
    friend TlSparseGeneralMatrix_EigenFloat operator*(const TlSparseSymmetricMatrix_EigenFloat& sm1,
                                                      const TlSparseSymmetricMatrix_EigenFloat& sm2);

    friend TlDenseVector_EigenFloat operator*(const TlSparseSymmetricMatrix_EigenFloat& mat,
                                         const TlDenseVector_EigenFloat& vtr);
};

#endif  // TL_SPARSE_SYMMETRIC_MATRIX_EIGEN_FLOAT_H
