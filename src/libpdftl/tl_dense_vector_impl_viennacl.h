#ifndef TL_DENSE_VECTOR_IMPL_VIENNACL_H
#define TL_DENSE_VECTOR_IMPL_VIENNACL_H

// IMPORTANT: Must be set prior to any ViennaCL includes if you want to use
// ViennaCL algorithms on Eigen objects
#define VIENNACL_WITH_EIGEN 1

#include "tl_dense_vector_impl_object.h"
#include "viennacl/vector.hpp"

class TlDenseGeneralMatrix_ImplViennaCL;

class TlDenseVector_ImplViennaCL : public TlDenseVector_ImplObject {
 public:
  typedef viennacl::vector<double> VectorDataType;

  // ---------------------------------------------------------------------------
  // constructor & destructor
  // ---------------------------------------------------------------------------
 public:
  TlDenseVector_ImplViennaCL(TlDenseVectorObject::index_type dim = 0);
  TlDenseVector_ImplViennaCL(const TlDenseVector_ImplViennaCL& rhs);
  virtual ~TlDenseVector_ImplViennaCL();

  // ---------------------------------------------------------------------------
  // properties
  // ---------------------------------------------------------------------------
 public:
  virtual TlDenseVectorObject::size_type getSize() const;
  virtual void resize(TlDenseVectorObject::index_type newSize);

 public:
  virtual double get(const TlDenseVectorObject::index_type i) const;
  virtual void set(const TlDenseVectorObject::index_type i, const double value);
  virtual void add(const TlDenseVectorObject::index_type i, const double value);
  virtual void mul(const TlDenseVectorObject::index_type i, const double value);

  // ---------------------------------------------------------------------------
  // operators
  // ---------------------------------------------------------------------------
 public:
  TlDenseVector_ImplViennaCL& operator=(const TlDenseVector_ImplViennaCL& rhs);

  TlDenseVector_ImplViennaCL& operator+=(const TlDenseVector_ImplViennaCL& rhs);
  TlDenseVector_ImplViennaCL& operator-=(const TlDenseVector_ImplViennaCL& rhs);
  TlDenseVector_ImplViennaCL& operator*=(const double rhs);
  TlDenseVector_ImplViennaCL& operator/=(const double rhs);

  // ---------------------------------------------------------------------------
  // operations
  // ---------------------------------------------------------------------------
 public:
  double sum() const;
  void sortByGreater();
  TlDenseVector_ImplViennaCL& dotInPlace(const TlDenseVector_ImplViennaCL& rhs);

  // ---------------------------------------------------------------------------
  // variables
  // ---------------------------------------------------------------------------
 protected:
  VectorDataType vector_;

  // ---------------------------------------------------------------------------
  // others
  // ---------------------------------------------------------------------------
  friend class TlDenseSymmetricMatrix_ImplViennaCL;

  friend TlDenseVector_ImplViennaCL operator+(
      const TlDenseVector_ImplViennaCL& rhs1,
      const TlDenseVector_ImplViennaCL& rhs2);
  friend TlDenseVector_ImplViennaCL operator-(
      const TlDenseVector_ImplViennaCL& rhs1,
      const TlDenseVector_ImplViennaCL& rhs2);

  friend TlDenseVector_ImplViennaCL operator*(
      const TlDenseGeneralMatrix_ImplViennaCL& mat,
      const TlDenseVector_ImplViennaCL& vec);
  friend TlDenseVector_ImplViennaCL operator*(
      const TlDenseVector_ImplViennaCL& rhs1, const double rhs2);
  friend TlDenseVector_ImplViennaCL operator*(
      const double rhs1, const TlDenseVector_ImplViennaCL& rhs2);

  friend TlDenseVector_ImplViennaCL operator*(
      const TlDenseGeneralMatrix_ImplViennaCL& mat,
      const TlDenseVector_ImplViennaCL& vec);
  friend TlDenseVector_ImplViennaCL operator*(
      const TlDenseVector_ImplViennaCL& vec,
      const TlDenseGeneralMatrix_ImplViennaCL& mat);
};

#endif  // TL_DENSE_VECTOR_IMPL_VIENNACL_H
