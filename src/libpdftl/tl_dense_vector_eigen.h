#ifndef TL_DENSE_VECTOR_EIGEN_H
#define TL_DENSE_VECTOR_EIGEN_H

#include "tl_dense_vector_object.h"
class TlDenseGeneralMatrix_Eigen;
class TlDenseSymmetricMatrix_Eigen;

class TlDenseVector_Eigen : public TlDenseVectorObject {
 public:
  explicit TlDenseVector_Eigen(TlDenseVectorObject::index_type size = 0);
  TlDenseVector_Eigen(const double* p,
                      const TlDenseVectorObject::size_type length);
  TlDenseVector_Eigen(const TlDenseVector_Eigen& rhs);

  virtual ~TlDenseVector_Eigen();

  // ---------------------------------------------------------------------------
  // operators
  // ---------------------------------------------------------------------------
 public:
  TlDenseVector_Eigen& operator=(const TlDenseVector_Eigen& rhs);

  TlDenseVector_Eigen& operator+=(const TlDenseVector_Eigen& rhs);
  TlDenseVector_Eigen& operator-=(const TlDenseVector_Eigen& rhs);
  TlDenseVector_Eigen& operator*=(const double rhs);
  TlDenseVector_Eigen& operator/=(const double rhs);

  // ---------------------------------------------------------------------------
  // operations
  // ---------------------------------------------------------------------------
  TlDenseVector_Eigen& dotInPlace(const TlDenseVector_Eigen& rhs);

  // ---------------------------------------------------------------------------
  // others
  // ---------------------------------------------------------------------------
  friend class TlDenseGeneralMatrix_Eigen;
  friend class TlDenseSymmetricMatrix_Eigen;

  friend TlDenseVector_Eigen operator+(const TlDenseVector_Eigen& rhs1,
                                       const TlDenseVector_Eigen& rhs2);
  friend TlDenseVector_Eigen operator-(const TlDenseVector_Eigen& rhs1,
                                       const TlDenseVector_Eigen& rhs2);
  friend TlDenseVector_Eigen operator*(const TlDenseVector_Eigen& rhs1,
                                       const double rhs2);
  friend TlDenseVector_Eigen operator*(const double rhs1,
                                       const TlDenseVector_Eigen& rhs2);

  friend TlDenseVector_Eigen operator*(
      const TlDenseGeneralMatrix_Eigen& rhs1,
      const TlDenseVector_Eigen& rhs2);
  friend TlDenseVector_Eigen operator*(
      const TlDenseVector_Eigen& rhs1,
      const TlDenseGeneralMatrix_Eigen& rhs2);
};

#endif  // TL_DENSE_VECTOR_EIGEN_H
