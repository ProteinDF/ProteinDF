#ifndef TL_DENSE_SYMMETRIC_MATRIX_EIGEN_H
#define TL_DENSE_SYMMETRIC_MATRIX_EIGEN_H

#include "tl_dense_symmetric_matrix_object.h"

class TlDenseGeneralMatrix_Eigen;
class TlDenseVector_Eigen;

class TlDenseSymmetricMatrix_Eigen : public TlDenseSymmetricMatrixObject {
  // ---------------------------------------------------------------------------
  // constructor & destructor
  // ---------------------------------------------------------------------------
 public:
  explicit TlDenseSymmetricMatrix_Eigen(
      const TlMatrixObject::index_type dim = 1);
  TlDenseSymmetricMatrix_Eigen(const TlDenseSymmetricMatrix_Eigen& rhs);
  TlDenseSymmetricMatrix_Eigen(const TlDenseGeneralMatrix_Eigen& rhs);
  virtual ~TlDenseSymmetricMatrix_Eigen();

  // ---------------------------------------------------------------------------
  // operators
  // ---------------------------------------------------------------------------
 public:
  TlDenseSymmetricMatrix_Eigen& operator=(
      const TlDenseSymmetricMatrix_Eigen& rhs);

  const TlDenseSymmetricMatrix_Eigen operator+(
      const TlDenseSymmetricMatrix_Eigen& rhs) const;
  const TlDenseSymmetricMatrix_Eigen operator-(
      const TlDenseSymmetricMatrix_Eigen& rhs) const;
  const TlDenseSymmetricMatrix_Eigen operator*(
      const TlDenseSymmetricMatrix_Eigen& rhs) const;

  TlDenseSymmetricMatrix_Eigen& operator+=(
      const TlDenseSymmetricMatrix_Eigen& rhs);
  TlDenseSymmetricMatrix_Eigen& operator-=(
      const TlDenseSymmetricMatrix_Eigen& rhs);
  TlDenseSymmetricMatrix_Eigen& operator*=(const double coef);
  TlDenseSymmetricMatrix_Eigen& operator/=(const double coef);
  TlDenseSymmetricMatrix_Eigen& operator*=(
      const TlDenseSymmetricMatrix_Eigen& rhs);

  // ---------------------------------------------------------------------------
  // operations
  // ---------------------------------------------------------------------------
  // double sum() const;
  const TlDenseSymmetricMatrix_Eigen& dotInPlace(
      const TlDenseSymmetricMatrix_Eigen& rhs);

  bool eig(TlDenseVector_Eigen* pEigVal,
           TlDenseGeneralMatrix_Eigen* pEigVec) const;
  bool diagonal(TlDenseVector_Eigen* pEigVal,
                TlDenseGeneralMatrix_Eigen* pEigVec) const {
    return this->eig(pEigVal, pEigVec);
  }

  TlDenseSymmetricMatrix_Eigen inverse() const;

  // ---------------------------------------------------------------------------
  // others
  // ---------------------------------------------------------------------------
  friend class TlDenseGeneralMatrix_Eigen;
  friend class TlDenseSymmetricMatrix_ViennaCL;
};

#endif  // TL_DENSE_SYMMETRIC_MATRIX_EIGEN_H
