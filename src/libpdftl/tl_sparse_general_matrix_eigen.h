#ifndef TL_SPARSE_GENERAL_MATRIX_EIGEN_H
#define TL_SPARSE_GENERAL_MATRIX_EIGEN_H

#include "tl_sparse_general_matrix_impl_eigen.h"
#include "tl_sparse_general_matrix_object.h"

class TlDenseGeneralMatrix_Eigen;
class TlDenseVector_Eigen;
class TlSparseSymmetricMatrix_Eigen;
class TlSparseSymmetricMatrix_ImplEigen;
class TlSparseGeneralMatrix_ViennaCL;

class TlSparseGeneralMatrix_Eigen : public TlSparseGeneralMatrixObject {
  // ---------------------------------------------------------------------------
  // constructor & destructor
  // ---------------------------------------------------------------------------
 public:
  explicit TlSparseGeneralMatrix_Eigen(
      const TlMatrixObject::index_type row = 1,
      const TlMatrixObject::index_type col = 1);
  TlSparseGeneralMatrix_Eigen(const TlSparseGeneralMatrix_Eigen& rhs);
  TlSparseGeneralMatrix_Eigen(const TlSparseSymmetricMatrix_Eigen& rhs);
  TlSparseGeneralMatrix_Eigen(const TlSparseGeneralMatrix_ImplEigen& rhs);
  TlSparseGeneralMatrix_Eigen(const TlDenseGeneralMatrix_Eigen& rhs);

  #ifdef HAVE_VIENNACL
  TlSparseGeneralMatrix_Eigen(const TlSparseGeneralMatrix_ViennaCL& rhs);
  #endif // HAVE_VIENNACL

  virtual ~TlSparseGeneralMatrix_Eigen();

  // ---------------------------------------------------------------------------
  // operators
  // ---------------------------------------------------------------------------
 public:
  TlSparseGeneralMatrix_Eigen& operator=(
      const TlSparseGeneralMatrix_Eigen& rhs);
  TlSparseGeneralMatrix_Eigen& operator+=(
      const TlSparseGeneralMatrix_Eigen& rhs);
  TlSparseGeneralMatrix_Eigen& operator-=(
      const TlSparseGeneralMatrix_Eigen& rhs);
  TlSparseGeneralMatrix_Eigen& operator*=(const double coef);

  // ---------------------------------------------------------------------------
  // others
  // ---------------------------------------------------------------------------
  friend class TlDenseGeneralMatrix_Eigen;
  friend class TlSparseSymmetricMatrix_Eigen;
  friend class TlSparseSymmetricMatrix_ImplEigen;
  friend class TlSparseGeneralMatrix_ViennaCL;

  // DM(G) = DM(G) * SM(G)
  friend TlDenseGeneralMatrix_Eigen operator*(
      const TlDenseGeneralMatrix_Eigen& mat1,
      const TlSparseGeneralMatrix_Eigen& mat2);
  // DM(G) = SM(G) * DM(G)
  friend TlDenseGeneralMatrix_Eigen operator*(
      const TlSparseGeneralMatrix_Eigen& mat1,
      const TlDenseGeneralMatrix_Eigen& mat2);

  // SM(G) = SM(G) * SM(S)
  friend TlSparseGeneralMatrix_Eigen operator*(
      const TlSparseGeneralMatrix_Eigen& sm1,
      const TlSparseSymmetricMatrix_Eigen& sm2);
  // SM(G) = SM(S) * SM(G)
  friend TlSparseGeneralMatrix_Eigen operator*(
      const TlSparseSymmetricMatrix_Eigen& sm1,
      const TlSparseGeneralMatrix_Eigen& sm2);
  // SM(G) = SM(G) * SM(G)
  friend TlSparseGeneralMatrix_Eigen operator*(
      const TlSparseGeneralMatrix_Eigen& sm1,
      const TlSparseGeneralMatrix_Eigen& sm2);

  // DV = SM(G) * DV
  friend TlDenseVector_Eigen operator*(const TlSparseGeneralMatrix_Eigen& mat,
                                       const TlDenseVector_Eigen& vtr);
};
#endif  // TL_SPARSE_GENERAL_MATRIX_EIGEN_H
