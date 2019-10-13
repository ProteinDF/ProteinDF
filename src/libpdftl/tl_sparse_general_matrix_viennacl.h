#ifndef TL_SPARSE_GENERAL_MATRIX_VIENNACL_H
#define TL_SPARSE_GENERAL_MATRIX_VIENNACL_H

#include "tl_sparse_general_matrix_impl_viennacl.h"
#include "tl_sparse_general_matrix_object.h"

class TlSparseSymmetricMatrix_ViennaCL;
class TlSparseSymmetricMatrix_ImplViennaCL;
class TlDenseGeneralMatrix_ViennaCL;
class TlDenseVector_ViennaCL;
class TlSparseGeneralMatrix_Eigen;

class TlSparseGeneralMatrix_ViennaCL : public TlSparseGeneralMatrixObject {
    // ---------------------------------------------------------------------------
    // constructor & destructor
    // ---------------------------------------------------------------------------
   public:
    explicit TlSparseGeneralMatrix_ViennaCL(
        const TlMatrixObject::index_type row = 1,
        const TlMatrixObject::index_type col = 1);
    TlSparseGeneralMatrix_ViennaCL(const TlSparseGeneralMatrix_ViennaCL& rhs);
    TlSparseGeneralMatrix_ViennaCL(const TlSparseSymmetricMatrix_ViennaCL& rhs);
    TlSparseGeneralMatrix_ViennaCL(
        const TlSparseGeneralMatrix_ImplViennaCL& rhs);
    TlSparseGeneralMatrix_ViennaCL(const TlDenseGeneralMatrix_ViennaCL& rhs);

#ifdef HAVE_EIGEN
    TlSparseGeneralMatrix_ViennaCL(const TlSparseGeneralMatrix_Eigen& rhs);
#endif  // HAVE_EIGEN

    virtual ~TlSparseGeneralMatrix_ViennaCL();

    // ---------------------------------------------------------------------------
    // operators
    // ---------------------------------------------------------------------------
   public:
    TlSparseGeneralMatrix_ViennaCL& operator=(
        const TlSparseGeneralMatrix_ViennaCL& rhs);
    TlSparseGeneralMatrix_ViennaCL& operator+=(
        const TlSparseGeneralMatrix_ViennaCL& rhs);
    TlSparseGeneralMatrix_ViennaCL& operator-=(
        const TlSparseGeneralMatrix_ViennaCL& rhs);
    TlSparseGeneralMatrix_ViennaCL& operator*=(const double coef);

    // ---------------------------------------------------------------------------
    // others
    // ---------------------------------------------------------------------------
    friend class TlSparseSymmetricMatrix_ViennaCL;
    friend class TlSparseSymmetricMatrix_ImplViennaCL;
    friend class TlSparseGeneralMatrix_Eigen;
    friend class TlDenseGeneralMatrix_ViennaCL;

    // SM(G) = SM(G) * SM(G)
    friend TlSparseGeneralMatrix_ViennaCL operator*(
        const TlSparseGeneralMatrix_ViennaCL& sm1,
        const TlSparseGeneralMatrix_ViennaCL& sm2);
    // SM(G) = SM(G) * SM(S)
    friend TlSparseGeneralMatrix_ViennaCL operator*(
        const TlSparseGeneralMatrix_ViennaCL& sm1,
        const TlSparseSymmetricMatrix_ViennaCL& sm2);
    // SM(G) = SM(S) * SM(G)
    friend TlSparseGeneralMatrix_ViennaCL operator*(
        const TlSparseSymmetricMatrix_ViennaCL& sm1,
        const TlSparseGeneralMatrix_ViennaCL& sm2);
    // DV = SM(G) * DV
    friend TlDenseVector_ViennaCL operator*(
        const TlSparseGeneralMatrix_ViennaCL& mat,
        const TlDenseVector_ViennaCL& vtr);
};
#endif  // TL_SPARSE_GENERAL_MATRIX_VIENNACL_H
