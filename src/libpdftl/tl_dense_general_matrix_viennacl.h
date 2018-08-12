#ifndef TL_DENSE_GENERAL_MATRIX_VIENNACL_H
#define TL_DENSE_GENERAL_MATRIX_VIENNACL_H

#include "tl_dense_general_matrix_object.h"

class TlDenseSymmetricMatrix_ViennaCL;
class TlDenseGeneralMatrix_ImplViennaCL;
class TlDenseVector_ViennaCL;

class TlDenseGeneralMatrix_Eigen;

class TlDenseGeneralMatrix_ViennaCL : public TlDenseGeneralMatrixObject {
  // ---------------------------------------------------------------------------
  // constructor & destructor
  // ---------------------------------------------------------------------------
 public:
  explicit TlDenseGeneralMatrix_ViennaCL(
      const TlMatrixObject::index_type row = 1,
      const TlMatrixObject::index_type col = 1);
  TlDenseGeneralMatrix_ViennaCL(const TlDenseGeneralMatrix_ViennaCL& rhs);
  TlDenseGeneralMatrix_ViennaCL(const TlDenseSymmetricMatrix_ViennaCL& rhs);
  TlDenseGeneralMatrix_ViennaCL(const TlDenseGeneralMatrix_ImplViennaCL& rhs);
  TlDenseGeneralMatrix_ViennaCL(const TlDenseGeneralMatrix_Eigen& rhs);
  virtual ~TlDenseGeneralMatrix_ViennaCL();

  // ---------------------------------------------------------------------------
  // operators
  // ---------------------------------------------------------------------------
 public:
  TlDenseGeneralMatrix_ViennaCL& operator=(
      const TlDenseGeneralMatrix_ViennaCL& rhs);

  const TlDenseGeneralMatrix_ViennaCL operator+(
      const TlDenseGeneralMatrix_ViennaCL& rhs) const;
  const TlDenseGeneralMatrix_ViennaCL operator-(
      const TlDenseGeneralMatrix_ViennaCL& rhs) const;
  const TlDenseGeneralMatrix_ViennaCL operator*(
      const TlDenseGeneralMatrix_ViennaCL& rhs) const;

  TlDenseGeneralMatrix_ViennaCL& operator+=(
      const TlDenseGeneralMatrix_ViennaCL& rhs);
  TlDenseGeneralMatrix_ViennaCL& operator-=(
      const TlDenseGeneralMatrix_ViennaCL& rhs);
  TlDenseGeneralMatrix_ViennaCL& operator*=(const double coef);
  TlDenseGeneralMatrix_ViennaCL& operator/=(const double coef);
  TlDenseGeneralMatrix_ViennaCL& operator*=(
      const TlDenseGeneralMatrix_ViennaCL& rhs);

  // ---------------------------------------------------------------------------
  // operations
  // ---------------------------------------------------------------------------
 public:
  virtual double sum() const;
  // virtual double trace() const;
  virtual double getRMS() const;
  // virtual double getMaxAbsoluteElement(
  //  TlMatrixObject::index_type* outRow,
  //  TlMatrixObject::index_type* outCol) const;

  const TlDenseGeneralMatrix_ViennaCL& dotInPlace(
      const TlDenseGeneralMatrix_ViennaCL& rhs);
  TlDenseGeneralMatrix_ViennaCL transpose() const;
  TlDenseGeneralMatrix_ViennaCL inverse() const;

  // ---------------------------------------------------------------------------
  // others
  // ---------------------------------------------------------------------------
  friend class TlDenseSymmetricMatrix_ViennaCL;

  friend TlDenseVector_ViennaCL operator*(
      const TlDenseGeneralMatrix_ViennaCL& rhs1,
      const TlDenseVector_ViennaCL& rhs2);
  friend TlDenseVector_ViennaCL operator*(
      const TlDenseVector_ViennaCL& rhs1,
      const TlDenseGeneralMatrix_ViennaCL& rhs2);
};

#endif  // TL_DENSE_GENERAL_MATRIX_VIENNACL_H
