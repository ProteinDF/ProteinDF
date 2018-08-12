#ifndef TL_DENSE_SYMMETRIC_MATRIX_VIENNACL_H
#define TL_DENSE_SYMMETRIC_MATRIX_VIENNACL_H

#include "tl_dense_symmetric_matrix_object.h"

class TlDenseGeneralMatrix_ViennaCL;
class TlDenseVector_ViennaCL;
class TlDenseSymmetricMatrix_Eigen;

class TlDenseSymmetricMatrix_ViennaCL : public TlDenseSymmetricMatrixObject {
  // ---------------------------------------------------------------------------
  // constructor & destructor
  // ---------------------------------------------------------------------------
 public:
  explicit TlDenseSymmetricMatrix_ViennaCL(
      const TlMatrixObject::index_type dim = 1);
  TlDenseSymmetricMatrix_ViennaCL(const TlDenseSymmetricMatrix_ViennaCL& rhs);
  TlDenseSymmetricMatrix_ViennaCL(const TlDenseGeneralMatrix_ViennaCL& rhs);
  TlDenseSymmetricMatrix_ViennaCL(const TlDenseSymmetricMatrix_Eigen& rhs);
  virtual ~TlDenseSymmetricMatrix_ViennaCL();

  // ---------------------------------------------------------------------------
  // operators
  // ---------------------------------------------------------------------------
 public:
  TlDenseSymmetricMatrix_ViennaCL& operator=(
      const TlDenseSymmetricMatrix_ViennaCL& rhs);
  TlDenseSymmetricMatrix_ViennaCL& operator=(
      const TlDenseSymmetricMatrix_Eigen& rhs);

  const TlDenseSymmetricMatrix_ViennaCL operator+(
      const TlDenseSymmetricMatrix_ViennaCL& rhs) const;
  const TlDenseSymmetricMatrix_ViennaCL operator-(
      const TlDenseSymmetricMatrix_ViennaCL& rhs) const;
  const TlDenseSymmetricMatrix_ViennaCL operator*(
      const TlDenseSymmetricMatrix_ViennaCL& rhs) const;

  TlDenseSymmetricMatrix_ViennaCL& operator+=(
      const TlDenseSymmetricMatrix_ViennaCL& rhs);
  TlDenseSymmetricMatrix_ViennaCL& operator-=(
      const TlDenseSymmetricMatrix_ViennaCL& rhs);
  TlDenseSymmetricMatrix_ViennaCL& operator*=(const double coef);
  TlDenseSymmetricMatrix_ViennaCL& operator/=(const double coef);
  TlDenseSymmetricMatrix_ViennaCL& operator*=(
      const TlDenseSymmetricMatrix_ViennaCL& rhs);

  // ---------------------------------------------------------------------------
  // operations
  // ---------------------------------------------------------------------------
  const TlDenseSymmetricMatrix_ViennaCL& dotInPlace(
      const TlDenseSymmetricMatrix_ViennaCL& rhs);

  bool eig(TlDenseVector_ViennaCL* pEigVal,
           TlDenseGeneralMatrix_ViennaCL* pEigVec) const;
  bool eig_QR(TlDenseVector_ViennaCL* pEigVal,
              TlDenseGeneralMatrix_ViennaCL* pEigVec) const;

  TlDenseSymmetricMatrix_ViennaCL inverse() const;

  // ---------------------------------------------------------------------------
  // Others
  // ---------------------------------------------------------------------------
  friend class TlDenseGeneralMatrix_ViennaCL;
};

#endif  // TL_DENSE_SYMMETRIC_MATRIX_VIENNACL_H
