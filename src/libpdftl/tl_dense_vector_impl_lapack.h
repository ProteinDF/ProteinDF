#ifndef TL_DENSE_VECTOR_IMPL_LAPACK_H
#define TL_DENSE_VECTOR_IMPL_LAPACK_H

#include <vector>
#include "tl_dense_vector_impl_object.h"

class TlDenseGeneralMatrix_ImplLapack;
class TlDenseSymmetricMatrix_ImplLapack;

class TlDenseVector_ImplLapack : public TlDenseVector_ImplObject {
  // ---------------------------------------------------------------------------
  // constructor & destructor
  // ---------------------------------------------------------------------------
 public:
  explicit TlDenseVector_ImplLapack(
      const TlDenseVectorObject::index_type size = 0);
  TlDenseVector_ImplLapack(const TlDenseVector_ImplLapack& rhs);
  TlDenseVector_ImplLapack(const std::vector<double>& rhs);
  virtual ~TlDenseVector_ImplLapack();

  // ---------------------------------------------------------------------------
  // properties
  // ---------------------------------------------------------------------------
 public:
  virtual TlDenseVectorObject::size_type getSize() const;
  virtual void resize(const TlDenseVectorObject::index_type newSize);

  virtual double get(const TlDenseVectorObject::index_type i) const;
  virtual void set(const TlDenseVectorObject::index_type i, const double value);
  virtual void add(const TlDenseVectorObject::index_type i, const double value);
  virtual void mul(const TlDenseVectorObject::index_type i, const double value);

  // ---------------------------------------------------------------------------
  // operators
  // ---------------------------------------------------------------------------
 public:
  TlDenseVector_ImplLapack& operator=(const TlDenseVector_ImplLapack& rhs);

  TlDenseVector_ImplLapack& operator+=(const TlDenseVector_ImplLapack& rhs);
  TlDenseVector_ImplLapack& operator-=(const TlDenseVector_ImplLapack& rhs);
  TlDenseVector_ImplLapack& operator*=(const double rhs);
  TlDenseVector_ImplLapack& operator/=(const double rhs);

  double operator*(const TlDenseVector_ImplLapack& rhs) const;
  // ---------------------------------------------------------------------------
  // operations
  // ---------------------------------------------------------------------------
 public:
  double* data();
  const double* data() const;

  virtual double sum() const;
  virtual void sortByGreater();

  TlDenseVector_ImplLapack& dotInPlace(const TlDenseVector_ImplLapack& rhs);

  // ---------------------------------------------------------------------------
  // protected
  // ---------------------------------------------------------------------------
 protected:
  void initialize(bool clearIfNeeded = true);

  // ---------------------------------------------------------------------------
  // variables
  // ---------------------------------------------------------------------------
 protected:
  TlDenseVectorObject::index_type size_;
  double* vector_;

  // ---------------------------------------------------------------------------
  // friend
  // ---------------------------------------------------------------------------
  friend class TlCommunicate;
  friend class TlDenseGeneralMatrix_ImplLapack;
  friend class TlDenseSymmetricMatrix_ImplLapack;

  friend TlDenseVector_ImplLapack operator+(
      const TlDenseVector_ImplLapack& rhs1,
      const TlDenseVector_ImplLapack& rhs2);
  friend TlDenseVector_ImplLapack operator-(
      const TlDenseVector_ImplLapack& rhs1,
      const TlDenseVector_ImplLapack& rhs2);

  friend TlDenseVector_ImplLapack operator*(
      const TlDenseGeneralMatrix_ImplLapack& mat,
      const TlDenseVector_ImplLapack& vec);
  friend TlDenseVector_ImplLapack operator*(
      const TlDenseVector_ImplLapack& vec,
      const TlDenseGeneralMatrix_ImplLapack& mat);
  friend TlDenseVector_ImplLapack operator*(
      const TlDenseSymmetricMatrix_ImplLapack& mat,
      const TlDenseVector_ImplLapack& vec);

  friend TlDenseVector_ImplLapack operator*(
      const TlDenseVector_ImplLapack& rhs1, const double rhs2);
  friend TlDenseVector_ImplLapack operator*(
      const double rhs1, const TlDenseVector_ImplLapack& rhs2);
};

#endif  // TL_DENSE_VECTOR_IMPL_LAPACK_H
