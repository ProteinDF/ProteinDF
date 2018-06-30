#ifndef TL_DENSE_VECTOR_VIENNACL_H
#define TL_DENSE_VECTOR_VIENNACL_H

#include "tl_dense_vector_object.h"

class TlDenseGeneralMatrix_ViennaCL;
class TlDenseSymmetricMatrix_ViennaCL;

class TlDenseVector_ViennaCL : public TlDenseVectorObject {
  // ---------------------------------------------------------------------------
  // constructor & destructor
  // ---------------------------------------------------------------------------
 public:
  explicit TlDenseVector_ViennaCL(
      const TlDenseVectorObject::index_type size = 0);
  TlDenseVector_ViennaCL(const TlDenseVector_ViennaCL& rhs);
  virtual ~TlDenseVector_ViennaCL();

  // ---------------------------------------------------------------------------
  // operators
  // ---------------------------------------------------------------------------
 public:
  TlDenseVector_ViennaCL& operator=(const TlDenseVector_ViennaCL& rhs);

  TlDenseVector_ViennaCL& operator+=(const TlDenseVector_ViennaCL& rhs);
  TlDenseVector_ViennaCL& operator-=(const TlDenseVector_ViennaCL& rhs);
  TlDenseVector_ViennaCL& operator*=(const double rhs);
  TlDenseVector_ViennaCL& operator/=(const double rhs);

  // ---------------------------------------------------------------------------
  // operations
  // ---------------------------------------------------------------------------
  TlDenseVector_ViennaCL& dotInPlace(const TlDenseVector_ViennaCL& rhs);

  // ---------------------------------------------------------------------------
  // others
  // ---------------------------------------------------------------------------
  friend class TlDenseSymmetricMatrix_ViennaCL;

  friend TlDenseVector_ViennaCL operator+(const TlDenseVector_ViennaCL& rhs1,
                                          const TlDenseVector_ViennaCL& rhs2);
  friend TlDenseVector_ViennaCL operator-(const TlDenseVector_ViennaCL& rhs1,
                                          const TlDenseVector_ViennaCL& rhs2);
  friend TlDenseVector_ViennaCL operator*(const TlDenseVector_ViennaCL& rhs1,
                                          const double rhs2);
  friend TlDenseVector_ViennaCL operator*(const double rhs1,
                                          const TlDenseVector_ViennaCL& rhs2);

  friend TlDenseVector_ViennaCL operator*(
      const TlDenseGeneralMatrix_ViennaCL& rhs1,
      const TlDenseVector_ViennaCL& rhs2);
  friend TlDenseVector_ViennaCL operator*(
      const TlDenseVector_ViennaCL& rhs1,
      const TlDenseGeneralMatrix_ViennaCL& rhs2);
};

#endif  // TL_DENSE_VECTOR_VIENNACL_H
