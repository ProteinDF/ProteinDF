#ifndef TL_DENSE_VECTOR_LAPACK_H
#define TL_DENSE_VECTOR_LAPACK_H

#include "tl_dense_vector_object.h"

class TlDenseGeneralMatrix_Lapack;
class TlDenseSymmetricMatrix_Lapack;

class TlDenseVector_Lapack : public TlDenseVectorObject {
  // ---------------------------------------------------------------------------
  // constructor & destructor
  // ---------------------------------------------------------------------------
 public:
  explicit TlDenseVector_Lapack(const TlDenseVectorObject::index_type size = 0);
  TlDenseVector_Lapack(const TlDenseVector_Lapack& rhs);
  // TlDenseVector_Lapack(const double* p,
  //                      const TlDenseVectorObject::size_type length);
  virtual ~TlDenseVector_Lapack();

  // ---------------------------------------------------------------------------
  // operators
  // ---------------------------------------------------------------------------
 public:
  TlDenseVector_Lapack& operator=(const TlDenseVector_Lapack& rhs);

  TlDenseVector_Lapack& operator+=(const TlDenseVector_Lapack& rhs);
  TlDenseVector_Lapack& operator-=(const TlDenseVector_Lapack& rhs);
  TlDenseVector_Lapack& operator*=(const double rhs);
  TlDenseVector_Lapack& operator/=(const double rhs);

  // ---------------------------------------------------------------------------
  // operations
  // ---------------------------------------------------------------------------
  TlDenseVector_Lapack& dotInPlace(const TlDenseVector_Lapack& rhs);

  // ---------------------------------------------------------------------------
  // others
  // ---------------------------------------------------------------------------
  friend class TlDenseSymmetricMatrix_Lapack;

  friend TlDenseVector_Lapack operator+(const TlDenseVector_Lapack& rhs1,
                                        const TlDenseVector_Lapack& rhs2);
  friend TlDenseVector_Lapack operator-(const TlDenseVector_Lapack& rhs1,
                                        const TlDenseVector_Lapack& rhs2);
  friend TlDenseVector_Lapack operator*(const TlDenseVector_Lapack& rhs1,
                                        const double rhs2);
  friend TlDenseVector_Lapack operator*(const double rhs1,
                                        const TlDenseVector_Lapack& rhs2);

  friend TlDenseVector_Lapack operator*(
      const TlDenseGeneralMatrix_Lapack& rhs1,
      const TlDenseVector_Lapack& rhs2);
  friend TlDenseVector_Lapack operator*(
      const TlDenseVector_Lapack& rhs1,
      const TlDenseGeneralMatrix_Lapack& rhs2);
};

#endif  // TL_DENSE_VECTOR_LAPACK_H
