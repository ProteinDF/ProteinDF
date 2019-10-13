#ifndef TL_SPARSE_SYMMETRIC_MATRIX_IMPL_VIENNACL_H
#define TL_SPARSE_SYMMETRIC_MATRIX_IMPL_VIENNACL_H

#include "tl_sparse_general_matrix_impl_viennacl.h"

#define VIENNACL_WITH_EIGEN 1

class TlSparseSymmetricMatrix_ImplEigen;
class TlDenseSymmetricMatrix_ImplViennaCL;

class TlSparseSymmetricMatrix_ImplViennaCL
    : public TlSparseGeneralMatrix_ImplViennaCL {
    // ---------------------------------------------------------------------------
    // constructor & destructor
    // ---------------------------------------------------------------------------
   public:
    explicit TlSparseSymmetricMatrix_ImplViennaCL(
        const TlMatrixObject::index_type dim = 0);
    TlSparseSymmetricMatrix_ImplViennaCL(
        const TlSparseSymmetricMatrix_ImplViennaCL& rhs);
    TlSparseSymmetricMatrix_ImplViennaCL(
        const TlSparseGeneralMatrix_ImplViennaCL& rhs);
    TlSparseSymmetricMatrix_ImplViennaCL(
        const TlDenseSymmetricMatrix_ImplViennaCL& rhs);

#ifdef HAVE_EIGEN
    TlSparseSymmetricMatrix_ImplViennaCL(
        const TlSparseSymmetricMatrix_ImplEigen& rhs);
#endif  // HAVE_EIGEN

    virtual ~TlSparseSymmetricMatrix_ImplViennaCL();

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
    TlSparseSymmetricMatrix_ImplViennaCL& operator=(
        const TlSparseSymmetricMatrix_ImplViennaCL& rhs);
    TlSparseSymmetricMatrix_ImplViennaCL& operator+=(
        const TlSparseSymmetricMatrix_ImplViennaCL& rhs);
    TlSparseSymmetricMatrix_ImplViennaCL& operator-=(
        const TlSparseSymmetricMatrix_ImplViennaCL& rhs);
    TlSparseSymmetricMatrix_ImplViennaCL& operator*=(const double coef);

    // ---------------------------------------------------------------------------
    // protected:
    // ---------------------------------------------------------------------------

    // ---------------------------------------------------------------------------
    // others
    // ---------------------------------------------------------------------------
    friend class TlSparseSymmetricMatrix_ImplEigen;

    // SM(G) = SM(G) * SM(S)
    friend TlSparseGeneralMatrix_ImplViennaCL operator*(
        const TlSparseGeneralMatrix_ImplViennaCL& sm1,
        const TlSparseSymmetricMatrix_ImplViennaCL& sm2);
    // SM(G) = SM(S) * SM(G)
    friend TlSparseGeneralMatrix_ImplViennaCL operator*(
        const TlSparseSymmetricMatrix_ImplViennaCL& sm1,
        const TlSparseGeneralMatrix_ImplViennaCL& sm2);
    // SM(G) = SM(S) * SM(S)
    friend TlSparseGeneralMatrix_ImplViennaCL operator*(
        const TlSparseSymmetricMatrix_ImplViennaCL& sm1,
        const TlSparseSymmetricMatrix_ImplViennaCL& sm2);

    friend TlDenseVector_ImplViennaCL operator*(
        const TlSparseSymmetricMatrix_ImplViennaCL& mat,
        const TlDenseVector_ImplViennaCL& vtr);
};

#endif  // TL_SPARSE_SYMMETRIC_MATRIX_IMPL_VIENNACL_H
