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
      const TlMatrixObject::index_type dim = 1,
      double const* const pBuf = NULL);
  TlDenseSymmetricMatrix_Lapack(const TlDenseSymmetricMatrix_Lapack& rhs);
  TlDenseSymmetricMatrix_Lapack(const TlDenseGeneralMatrix_Lapack& rhs);

  virtual ~TlDenseSymmetricMatrix_Lapack();

  operator std::vector<double>() const;

  // ---------------------------------------------------------------------------
  // properties
  // ---------------------------------------------------------------------------
 public:
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
  // TlDenseSymmetricMatrix_Lapack& operator*=(
  //     const TlDenseSymmetricMatrix_Lapack& rhs);

  // ---------------------------------------------------------------------------
  // operations
  // ---------------------------------------------------------------------------
  const TlDenseSymmetricMatrix_Lapack& dotInPlace(
      const TlDenseSymmetricMatrix_Lapack& rhs);

  bool eig(TlDenseVector_Lapack* pEigVal,
           TlDenseGeneralMatrix_Lapack* pEigVec) const;

  TlDenseSymmetricMatrix_Lapack inverse() const;

 public:
  double* data();
  const double* data() const;

  // ---------------------------------------------------------------------------
  // I/O
  // ---------------------------------------------------------------------------
 public:
  void dump(TlDenseVector_Lapack* v) const;
  void restore(const TlDenseVector_Lapack& v);

  // ---------------------------------------------------------------------------
  // protected
  // ---------------------------------------------------------------------------
 protected:
  virtual TlMatrixObject::size_type getNumOfElements() const;

  // ---------------------------------------------------------------------------
  // friend
  // ---------------------------------------------------------------------------
  friend class TlCommunicate;
  friend class TlDenseGeneralMatrix_Lapack;

  // matrix x matrix
  friend TlDenseGeneralMatrix_Lapack operator*(
      const TlDenseSymmetricMatrix_Lapack& rhs1,
      const TlDenseGeneralMatrix_Lapack& rhs2);
  friend TlDenseGeneralMatrix_Lapack operator*(
      const TlDenseGeneralMatrix_Lapack& rhs1,
      const TlDenseSymmetricMatrix_Lapack& rhs2);

  // matrix x vector
  friend TlDenseVector_Lapack operator*(
      const TlDenseSymmetricMatrix_Lapack& rhs1,
      const TlDenseVector_Lapack& rhs2);
  friend TlDenseVector_Lapack operator*(
      const TlDenseVector_Lapack& rhs1,
      const TlDenseSymmetricMatrix_Lapack& rhs2);
};

// ---------------------------------------------------------------------------
// arithmetic
// ---------------------------------------------------------------------------
TlDenseSymmetricMatrix_Lapack operator*(
    const double coef, const TlDenseSymmetricMatrix_Lapack& matrix);
TlDenseSymmetricMatrix_Lapack operator*(
    const TlDenseSymmetricMatrix_Lapack& matrix, const double coef);

#endif  // TL_DENSE_SYMMETRIC_MATRIX_LAPACK_H
