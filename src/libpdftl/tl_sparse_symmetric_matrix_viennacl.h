#ifndef TL_SPARSE_SYMMETRIC_MATRIX_VIENNACL_H
#define TL_SPARSE_SYMMETRIC_MATRIX_VIENNACL_H

#include "tl_sparse_symmetric_matrix_object.h"

class TlSparseGeneralMatrix_ViennaCL;
class TlSparseGeneralMatrix_ImplViennaCL;
class TlSparseSymmetricMatrix_Eigen;
class TlDenseSymmetricMatrix_ViennaCL;

class TlSparseSymmetricMatrix_ViennaCL : public TlSparseSymmetricMatrixObject {
    // ---------------------------------------------------------------------------
    // constructor & destructor
    // ---------------------------------------------------------------------------
   public:
    explicit TlSparseSymmetricMatrix_ViennaCL(
        const TlMatrixObject::index_type dim = 1);
    TlSparseSymmetricMatrix_ViennaCL(
        const TlSparseSymmetricMatrix_ViennaCL& rhs);
    TlSparseSymmetricMatrix_ViennaCL(const TlSparseGeneralMatrix_ViennaCL& rhs);
    TlSparseSymmetricMatrix_ViennaCL(
        const TlDenseSymmetricMatrix_ViennaCL& rhs);

#ifdef HAVE_EIGEN
    TlSparseSymmetricMatrix_ViennaCL(const TlSparseSymmetricMatrix_Eigen& rhs);
#endif  // HAVE_EIGEN

    virtual ~TlSparseSymmetricMatrix_ViennaCL();

    // ---------------------------------------------------------------------------
    // operators
    // ---------------------------------------------------------------------------
   public:
    TlSparseSymmetricMatrix_ViennaCL& operator=(
        const TlSparseSymmetricMatrix_ViennaCL& rhs);
    TlSparseSymmetricMatrix_ViennaCL& operator+=(
        const TlSparseSymmetricMatrix_ViennaCL& rhs);
    TlSparseSymmetricMatrix_ViennaCL& operator-=(
        const TlSparseSymmetricMatrix_ViennaCL& rhs);
    TlSparseSymmetricMatrix_ViennaCL& operator*=(const double coef);

    // ---------------------------------------------------------------------------
    // others
    // ---------------------------------------------------------------------------
    friend class TlSparseGeneralMatrix_ViennaCL;
    friend class TlSparseGeneralMatrix_ImplViennaCL;
    friend class TlSparseSymmetricMatrix_Eigen;
    friend class TlDenseSymmetricMatrix_ViennaCL;

    friend TlSparseGeneralMatrix_ViennaCL operator*(
        const TlSparseGeneralMatrix_ViennaCL& sm1,
        const TlSparseSymmetricMatrix_ViennaCL& sm2);
    friend TlSparseGeneralMatrix_ViennaCL operator*(
        const TlSparseSymmetricMatrix_ViennaCL& sm1,
        const TlSparseGeneralMatrix_ViennaCL& sm2);
    friend TlSparseGeneralMatrix_ViennaCL operator*(
        const TlSparseSymmetricMatrix_ViennaCL& sm1,
        const TlSparseSymmetricMatrix_ViennaCL& sm2);
};

#endif  // TL_SPARSE_SYMMETRIC_MATRIX_VIENNACL_H
