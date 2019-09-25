#ifndef TL_SPARSE_SYMMETRIC_MATRIX_IMPL_EIGEN_H
#define TL_SPARSE_SYMMETRIC_MATRIX_IMPL_EIGEN_H

#include "tl_sparse_general_matrix_impl_eigen.h"

class TlDenseSymmetricMatrix_ImplEigen;
class TlSparseSymmetricMatrix_ImplViennaCL;

class TlSparseSymmetricMatrix_ImplEigen
    : public TlSparseGeneralMatrix_ImplEigen {
    // ---------------------------------------------------------------------------
    // constructor & destructor
    // ---------------------------------------------------------------------------
   public:
    explicit TlSparseSymmetricMatrix_ImplEigen(
        const TlMatrixObject::index_type dim = 0);
    TlSparseSymmetricMatrix_ImplEigen(
        const TlSparseSymmetricMatrix_ImplEigen& rhs);
    TlSparseSymmetricMatrix_ImplEigen(
        const TlSparseGeneralMatrix_ImplEigen& rhs);
    TlSparseSymmetricMatrix_ImplEigen(
        const TlDenseSymmetricMatrix_ImplEigen& rhs);

#ifdef HAVE_VIENNACL
    TlSparseSymmetricMatrix_ImplEigen(
        const TlSparseSymmetricMatrix_ImplViennaCL& rhs);
#endif  // HAVE_VIENNACL

    virtual ~TlSparseSymmetricMatrix_ImplEigen();

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
    // operators
    // ---------------------------------------------------------------------------
    TlSparseSymmetricMatrix_ImplEigen& operator+=(
        const TlSparseSymmetricMatrix_ImplEigen& sm);
    TlSparseSymmetricMatrix_ImplEigen& operator-=(
        const TlSparseSymmetricMatrix_ImplEigen& sm);

    // ---------------------------------------------------------------------------
    // others
    // ---------------------------------------------------------------------------
    friend class TlDenseSymmetricMatrix_ImplEigen;
    friend class TlSparseSymmetricMatrix_ImplViennaCL;

    // DM(G) = DM(G) * SM(S)
    friend TlDenseGeneralMatrix_ImplEigen operator*(
        const TlDenseGeneralMatrix_ImplEigen& mat1,
        const TlSparseSymmetricMatrix_ImplEigen& mat2);
    // DM(G) = SM(S) * DM(G)
    friend TlDenseGeneralMatrix_ImplEigen operator*(
        const TlSparseSymmetricMatrix_ImplEigen& mat1,
        const TlDenseGeneralMatrix_ImplEigen& mat2);

    // SM(G) = SM(G) * SM(S)
    friend TlSparseGeneralMatrix_ImplEigen operator*(
        const TlSparseGeneralMatrix_ImplEigen& mat1,
        const TlSparseSymmetricMatrix_ImplEigen& mat2);
    // SM(G) = SM(S) * SM(G)
    friend TlSparseGeneralMatrix_ImplEigen operator*(
        const TlSparseSymmetricMatrix_ImplEigen& mat1,
        const TlSparseGeneralMatrix_ImplEigen& mat2);
    // SM(G) = SM(S) * SM(S)
    friend TlSparseGeneralMatrix_ImplEigen operator*(
        const TlSparseSymmetricMatrix_ImplEigen& sm1,
        const TlSparseSymmetricMatrix_ImplEigen& sm2);

    friend TlDenseVector_ImplEigen operator*(
        const TlSparseSymmetricMatrix_ImplEigen& mat,
        const TlDenseVector_ImplEigen& vtr);

   public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW  // Eigen macro
};

#endif  // TL_SPARSE_SYMMETRIC_MATRIX_IMPL_EIGEN_H
