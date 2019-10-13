#ifndef TL_SPARSE_SYMMETRIC_MATRIX_EIGEN_H
#define TL_SPARSE_SYMMETRIC_MATRIX_EIGEN_H

#include "tl_sparse_symmetric_matrix_object.h"

class TlDenseGeneralMatrix_Eigen;
class TlDenseSymmetricMatrix_Eigen;
class TlDenseVector_Eigen;
class TlSparseGeneralMatrix_Eigen;
class TlSparseGeneralMatrix_ImplEigen;
class TlSparseSymmetricMatrix_ViennaCL;

class TlSparseSymmetricMatrix_Eigen : public TlSparseSymmetricMatrixObject {
    // ---------------------------------------------------------------------------
    // constructor & destructor
    // ---------------------------------------------------------------------------
   public:
    explicit TlSparseSymmetricMatrix_Eigen(
        const TlMatrixObject::index_type dim = 1);
    TlSparseSymmetricMatrix_Eigen(const TlSparseSymmetricMatrix_Eigen& rhs);
    TlSparseSymmetricMatrix_Eigen(const TlSparseGeneralMatrix_Eigen& rhs);
    TlSparseSymmetricMatrix_Eigen(const TlDenseSymmetricMatrix_Eigen& rhs);

#ifdef HAVE_VIENNACL
    TlSparseSymmetricMatrix_Eigen(const TlSparseSymmetricMatrix_ViennaCL& rhs);
#endif  // HAVE_VIENNACL

    virtual ~TlSparseSymmetricMatrix_Eigen();

    // ---------------------------------------------------------------------------
    // operators
    // ---------------------------------------------------------------------------
   public:
    TlSparseSymmetricMatrix_Eigen& operator=(
        const TlSparseSymmetricMatrix_Eigen& rhs);

    TlSparseSymmetricMatrix_Eigen& operator+=(
        const TlSparseSymmetricMatrix_Eigen& sm);
    TlSparseSymmetricMatrix_Eigen& operator-=(
        const TlSparseSymmetricMatrix_Eigen& sm);

    // ---------------------------------------------------------------------------
    // others
    // ---------------------------------------------------------------------------
    friend class TlDenseSymmetricMatrix_Eigen;
    friend class TlSparseGeneralMatrix_Eigen;
    friend class TlSparseGeneralMatrix_ImplEigen;
    friend class TlSparseSymmetricMatrix_ViennaCL;

    // DM(G) = DM(G) * SM(S)
    friend TlDenseGeneralMatrix_Eigen operator*(
        const TlDenseGeneralMatrix_Eigen& mat1,
        const TlSparseSymmetricMatrix_Eigen& mat2);
    // DM(G) = SM(S) * DM(G)
    friend TlDenseGeneralMatrix_Eigen operator*(
        const TlSparseSymmetricMatrix_Eigen& mat1,
        const TlDenseGeneralMatrix_Eigen& mat2);

    // SM(G) = SM(G) * SM(S)
    friend TlSparseGeneralMatrix_Eigen operator*(
        const TlSparseGeneralMatrix_Eigen& sm1,
        const TlSparseSymmetricMatrix_Eigen& sm2);
    // SM(G) = SM(S) * SM(G)
    friend TlSparseGeneralMatrix_Eigen operator*(
        const TlSparseSymmetricMatrix_Eigen& sm1,
        const TlSparseGeneralMatrix_Eigen& sm2);
    // SM(G) = SM(S) * SM(S)
    friend TlSparseGeneralMatrix_Eigen operator*(
        const TlSparseSymmetricMatrix_Eigen& sm1,
        const TlSparseSymmetricMatrix_Eigen& sm2);

    friend TlDenseVector_Eigen operator*(
        const TlSparseSymmetricMatrix_Eigen& mat,
        const TlDenseVector_Eigen& vtr);
};

#endif  // TL_SPARSE_SYMMETRIC_MATRIX_EIGEN_H
