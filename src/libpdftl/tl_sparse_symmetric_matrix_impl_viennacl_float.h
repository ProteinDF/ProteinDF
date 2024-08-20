#ifndef TL_SPARSE_SYMMETRIC_MATRIX_IMPL_VIENNACL_H
#define TL_SPARSE_SYMMETRIC_MATRIX_IMPL_VIENNACL_H

#include "tl_sparse_general_matrix_impl_viennacl_float.h"

#define VIENNACL_WITH_EIGEN 1

class TlSparseSymmetricMatrix_ImplEigenFloat;
class TlDenseSymmetricMatrix_ImplViennaCLFloat;

class TlSparseSymmetricMatrix_ImplViennaCLFloat : public TlSparseGeneralMatrix_ImplViennaCLFloat {
    // ---------------------------------------------------------------------------
    // constructor & destructor
    // ---------------------------------------------------------------------------
public:
    explicit TlSparseSymmetricMatrix_ImplViennaCLFloat(const TlMatrixObject::index_type dim = 0);
    TlSparseSymmetricMatrix_ImplViennaCLFloat(const TlSparseSymmetricMatrix_ImplViennaCLFloat& rhs);
    TlSparseSymmetricMatrix_ImplViennaCLFloat(const TlSparseGeneralMatrix_ImplViennaCLFloat& rhs);
    TlSparseSymmetricMatrix_ImplViennaCLFloat(const TlDenseSymmetricMatrix_ImplViennaCLFloat& rhs);

#ifdef HAVE_EIGEN
    TlSparseSymmetricMatrix_ImplViennaCLFloat(const TlSparseSymmetricMatrix_ImplEigenFloat& rhs);
#endif  // HAVE_EIGEN

    virtual ~TlSparseSymmetricMatrix_ImplViennaCLFloat();

    // ---------------------------------------------------------------------------
    // properties
    // ---------------------------------------------------------------------------
    virtual double get(const TlMatrixObject::index_type row,
                       const TlMatrixObject::index_type col) const;

    virtual void set(const TlMatrixObject::index_type row,
                     const TlMatrixObject::index_type col, const double value);

    virtual void add(const TlMatrixObject::index_type row,
                     const TlMatrixObject::index_type col, const double value);

    virtual void mul(const TlMatrixObject::index_type row,
                     const TlMatrixObject::index_type col, const double value);
    // ---------------------------------------------------------------------------
    // operator
    // ---------------------------------------------------------------------------
    TlSparseSymmetricMatrix_ImplViennaCLFloat& operator=(
        const TlSparseSymmetricMatrix_ImplViennaCLFloat& rhs);
    TlSparseSymmetricMatrix_ImplViennaCLFloat& operator+=(
        const TlSparseSymmetricMatrix_ImplViennaCLFloat& rhs);
    TlSparseSymmetricMatrix_ImplViennaCLFloat& operator-=(
        const TlSparseSymmetricMatrix_ImplViennaCLFloat& rhs);
    TlSparseSymmetricMatrix_ImplViennaCLFloat& operator*=(const float coef);

    // ---------------------------------------------------------------------------
    // protected:
    // ---------------------------------------------------------------------------

    // ---------------------------------------------------------------------------
    // others
    // ---------------------------------------------------------------------------
    friend class TlSparseSymmetricMatrix_ImplEigenFloat;

    // SM(G) = SM(G) * SM(S)
    friend TlSparseGeneralMatrix_ImplViennaCLFloat operator*(
        const TlSparseGeneralMatrix_ImplViennaCLFloat& sm1,
        const TlSparseSymmetricMatrix_ImplViennaCLFloat& sm2);
    // SM(G) = SM(S) * SM(G)
    friend TlSparseGeneralMatrix_ImplViennaCLFloat operator*(
        const TlSparseSymmetricMatrix_ImplViennaCLFloat& sm1,
        const TlSparseGeneralMatrix_ImplViennaCLFloat& sm2);
    // SM(G) = SM(S) * SM(S)
    friend TlSparseGeneralMatrix_ImplViennaCLFloat operator*(
        const TlSparseSymmetricMatrix_ImplViennaCLFloat& sm1,
        const TlSparseSymmetricMatrix_ImplViennaCLFloat& sm2);

    friend TlDenseVector_ImplViennaCLFloat operator*(
        const TlSparseSymmetricMatrix_ImplViennaCLFloat& mat,
        const TlDenseVector_ImplViennaCLFloat& vtr);
};

#endif  // TL_SPARSE_SYMMETRIC_MATRIX_IMPL_VIENNACL_H
