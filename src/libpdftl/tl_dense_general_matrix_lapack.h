#ifndef TL_DENSE_GENERAL_MATRIX_LAPACK_H
#define TL_DENSE_GENERAL_MATRIX_LAPACK_H

#include "tl_dense_general_matrix_impl_lapack.h"
#include "tl_dense_general_matrix_object.h"

class TlDenseSymmetricMatrix_Lapack;
class TlDenseVector_Lapack;

class TlDenseGeneralMatrix_Lapack : public TlDenseGeneralMatrixObject {
  // ---------------------------------------------------------------------------
  // constructor & destructor
  // ---------------------------------------------------------------------------
 public:
  explicit TlDenseGeneralMatrix_Lapack(
      const TlMatrixObject::index_type row = 1,
      const TlMatrixObject::index_type col = 1);
  TlDenseGeneralMatrix_Lapack(const TlDenseGeneralMatrix_Lapack& rhs);
  TlDenseGeneralMatrix_Lapack(const TlDenseSymmetricMatrix_Lapack& rhs);
  TlDenseGeneralMatrix_Lapack(const TlDenseGeneralMatrix_ImplLapack& rhs);

  virtual ~TlDenseGeneralMatrix_Lapack();

  // ---------------------------------------------------------------------------
  // operators
  // ---------------------------------------------------------------------------
 public:
  TlDenseGeneralMatrix_Lapack& operator=(
      const TlDenseGeneralMatrix_Lapack& rhs);

  const TlDenseGeneralMatrix_Lapack operator+(
      const TlDenseGeneralMatrix_Lapack& rhs) const;
  const TlDenseGeneralMatrix_Lapack operator-(
      const TlDenseGeneralMatrix_Lapack& rhs) const;
  const TlDenseGeneralMatrix_Lapack operator*(
      const TlDenseGeneralMatrix_Lapack& rhs) const;

  TlDenseGeneralMatrix_Lapack& operator+=(
      const TlDenseGeneralMatrix_Lapack& rhs);
  TlDenseGeneralMatrix_Lapack& operator-=(
      const TlDenseGeneralMatrix_Lapack& rhs);
  TlDenseGeneralMatrix_Lapack& operator*=(const double coef);
  TlDenseGeneralMatrix_Lapack& operator/=(const double coef);
  TlDenseGeneralMatrix_Lapack& operator*=(
      const TlDenseGeneralMatrix_Lapack& rhs);

  // ---------------------------------------------------------------------------
  // operations
  // ---------------------------------------------------------------------------
 public:
  double sum() const;
  double getRMS() const;
  double getMaxAbsoluteElement(TlMatrixObject::index_type* outRow,
                               TlMatrixObject::index_type* outCol) const;

  const TlDenseGeneralMatrix_Lapack& dotInPlace(
      const TlDenseGeneralMatrix_Lapack& rhs);
  TlDenseGeneralMatrix_Lapack transpose() const;
  TlDenseGeneralMatrix_Lapack inverse() const;

  // ---------------------------------------------------------------------------
  // others
  // ---------------------------------------------------------------------------
  friend class TlDenseSymmetricMatrix_Lapack;

  friend TlDenseVector_Lapack operator*(const TlDenseGeneralMatrix_Lapack& rhs1,
                                        const TlDenseVector_Lapack& rhs2);
  friend TlDenseVector_Lapack operator*(
      const TlDenseVector_Lapack& rhs1,
      const TlDenseGeneralMatrix_Lapack& rhs2);
};

#endif  // TL_DENSE_GENERAL_MATRIX_LAPACK_H
