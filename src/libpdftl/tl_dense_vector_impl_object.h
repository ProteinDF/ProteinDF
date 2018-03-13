#ifndef TL_DENSE_VECTOR_IMPL_OBJECT_H
#define TL_DENSE_VECTOR_IMPL_OBJECT_H

#include "TlLogging.h"
#include "tl_dense_vector_object.h"

class TlDenseVector_ImplObject {
  // ---------------------------------------------------------------------------
  // constructor & destructor
  // ---------------------------------------------------------------------------
 public:
  TlDenseVector_ImplObject();
  virtual ~TlDenseVector_ImplObject();

  // ---------------------------------------------------------------------------
  // properties
  // ---------------------------------------------------------------------------
 public:
  virtual TlDenseVectorObject::index_type getSize() const = 0;
  virtual void resize(const TlDenseVectorObject::index_type newSize) = 0;

  virtual double get(const TlDenseVectorObject::index_type i) const = 0;
  virtual void set(const TlDenseVectorObject::index_type i,
                   const double value) = 0;
  virtual void add(const TlDenseVectorObject::index_type i,
                   const double value) = 0;

  // ---------------------------------------------------------------------------
  // operations
  // ---------------------------------------------------------------------------
  virtual double sum() const;
  virtual double norm() const;
  virtual double norm2() const;
  virtual TlDenseVectorObject::index_type argmax(
      const TlDenseVectorObject::index_type& begin,
      const TlDenseVectorObject::index_type& end) const;
  virtual void sortByGreater() = 0;

  // ---------------------------------------------------------------------------
  // variables
  // ---------------------------------------------------------------------------
 protected:
  TlLogging& log_;
};

#endif  // TL_DENSE_VECTOR_IMPL_OBJECT_H
