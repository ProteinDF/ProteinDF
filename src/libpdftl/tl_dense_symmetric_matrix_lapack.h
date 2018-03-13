#ifndef TL_DENSE_SYMMETRIC_MATRIX_LAPACK_H
#define TL_DENSE_SYMMETRIC_MATRIX_LAPACK_H

#include "tl_dense_symmetric_matrix_object.h"

class TlDenseGeneralMatrix_Lapack;
class TlDenseVector_Lapack;

class TlDenseSymmetricMatrix_Lapack : public TlDenseSymmetricMatrixObject {
  // ---------------------------------------------------------------------------
  // constructor & destructor
  // ---------------------------------------------------------------------------
 public:
  explicit TlDenseSymmetricMatrix_Lapack(
      const TlMatrixObject::index_type dim = 1);
  TlDenseSymmetricMatrix_Lapack(const TlDenseSymmetricMatrix_Lapack& rhs);
  TlDenseSymmetricMatrix_Lapack(const TlDenseGeneralMatrix_Lapack& rhs);
  virtual ~TlDenseSymmetricMatrix_Lapack();

  // ---------------------------------------------------------------------------
  // operators
  // ---------------------------------------------------------------------------
 public:
  // operator TlDenseGeneralMatrix_Lapack() const;

  TlDenseSymmetricMatrix_Lapack& operator=(
      const TlDenseSymmetricMatrix_Lapack& rhs);

  const TlDenseSymmetricMatrix_Lapack operator+(
      const TlDenseSymmetricMatrix_Lapack& rhs) const;
  const TlDenseSymmetricMatrix_Lapack operator-(
      const TlDenseSymmetricMatrix_Lapack& rhs) const;
  const TlDenseGeneralMatrix_Lapack operator*(
      const TlDenseSymmetricMatrix_Lapack& rhs) const;

  TlDenseSymmetricMatrix_Lapack& operator+=(
      const TlDenseSymmetricMatrix_Lapack& rhs);
  TlDenseSymmetricMatrix_Lapack& operator-=(
      const TlDenseSymmetricMatrix_Lapack& rhs);
  TlDenseSymmetricMatrix_Lapack& operator*=(const double coef);
  TlDenseSymmetricMatrix_Lapack& operator/=(const double coef);
  TlDenseSymmetricMatrix_Lapack& operator*=(
      const TlDenseSymmetricMatrix_Lapack& rhs);

  // ---------------------------------------------------------------------------
  // operations
  // ---------------------------------------------------------------------------
  const TlDenseSymmetricMatrix_Lapack& dotInPlace(
      const TlDenseSymmetricMatrix_Lapack& rhs);

  bool eig(TlDenseVector_Lapack* pEigVal,
           TlDenseGeneralMatrix_Lapack* pEigVec) const;

  TlDenseSymmetricMatrix_Lapack inverse() const;

  // ---------------------------------------------------------------------------
  // others
  // ---------------------------------------------------------------------------
  friend class TlDenseGeneralMatrix_Lapack;
};

#endif  // TL_DENSE_SYMMETRIC_MATRIX_LAPACK_H
