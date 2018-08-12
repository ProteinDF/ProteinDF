#ifndef TL_DENSE_GENERAL_MATRIX_EIGEN_H
#define TL_DENSE_GENERAL_MATRIX_EIGEN_H

#include "tl_dense_general_matrix_object.h"

class TlDenseSymmetricMatrix_Eigen;
class TlDenseVector_Eigen;
class TlDenseGeneralMatrix_ImplEigen;
class TlDenseGeneralMatrix_ViennaCL;

class TlDenseGeneralMatrix_Eigen : public TlDenseGeneralMatrixObject {
  // ---------------------------------------------------------------------------
  // constructor & destructor
  // ---------------------------------------------------------------------------
 public:
  explicit TlDenseGeneralMatrix_Eigen(const TlMatrixObject::index_type row = 1,
                                      const TlMatrixObject::index_type col = 1);
  TlDenseGeneralMatrix_Eigen(const TlDenseGeneralMatrix_Eigen& rhs);
  TlDenseGeneralMatrix_Eigen(const TlDenseSymmetricMatrix_Eigen& rhs);
  TlDenseGeneralMatrix_Eigen(const TlDenseGeneralMatrix_ImplEigen& rhs);
  virtual ~TlDenseGeneralMatrix_Eigen();

  // ---------------------------------------------------------------------------
  // operators
  // ---------------------------------------------------------------------------
 public:
  TlDenseGeneralMatrix_Eigen& operator=(const TlDenseGeneralMatrix_Eigen& rhs);

  const TlDenseGeneralMatrix_Eigen operator+(
      const TlDenseGeneralMatrix_Eigen& rhs) const;
  const TlDenseGeneralMatrix_Eigen operator-(
      const TlDenseGeneralMatrix_Eigen& rhs) const;
  const TlDenseGeneralMatrix_Eigen operator*(
      const TlDenseGeneralMatrix_Eigen& rhs) const;

  TlDenseGeneralMatrix_Eigen& operator+=(const TlDenseGeneralMatrix_Eigen& rhs);
  TlDenseGeneralMatrix_Eigen& operator-=(const TlDenseGeneralMatrix_Eigen& rhs);
  TlDenseGeneralMatrix_Eigen& operator*=(const double coef);
  TlDenseGeneralMatrix_Eigen& operator/=(const double coef);
  TlDenseGeneralMatrix_Eigen& operator*=(const TlDenseGeneralMatrix_Eigen& rhs);

  // ---------------------------------------------------------------------------
  // operations
  // ---------------------------------------------------------------------------
 public:
  double sum() const;
  double getRMS() const;

  const TlDenseGeneralMatrix_Eigen& dotInPlace(
      const TlDenseGeneralMatrix_Eigen& rhs);
  TlDenseGeneralMatrix_Eigen transpose() const;
  TlDenseGeneralMatrix_Eigen inverse() const;

  // ---------------------------------------------------------------------------
  // others
  // ---------------------------------------------------------------------------
  friend class TlDenseSymmetricMatrix_Eigen;
  friend class TlDenseGeneralMatrix_ViennaCL;

  friend TlDenseVector_Eigen operator*(const TlDenseGeneralMatrix_Eigen& rhs1,
                                       const TlDenseVector_Eigen& rhs2);
  friend TlDenseVector_Eigen operator*(const TlDenseVector_Eigen& rhs1,
                                       const TlDenseGeneralMatrix_Eigen& rhs2);
};

#endif  // TL_DENSE_GENERAL_MATRIX_EIGEN_H
