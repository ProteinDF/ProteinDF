#ifndef TL_DENSE_VECTOR_IMPL_EIGEN_H
#define TL_DENSE_VECTOR_IMPL_EIGEN_H

#include <vector>
#include <Eigen/Core>

#include "tl_dense_vector_impl_object.h"

#if __cplusplus >= 201103L
#include <mutex>
#endif  // __cplusplus

class TlDenseGeneralMatrix_ImplEigen;
class TlDenseSymmetricMatrix_ImplEigen;
class TlSparseGeneralMatrix_ImplEigen;
class TlSparseSymmetricMatrix_ImplEigen;
class TlDenseVector_ImplViennaCL;

class TlDenseVector_ImplEigen : public TlDenseVector_ImplObject {
 public:
  typedef Eigen::VectorXd VectorDataType;  // Matrix<double, Dynamic, 1>
  typedef Eigen::Map<VectorDataType> MapType;
  typedef Eigen::Map<const VectorDataType> MapTypeConst;

  // ---------------------------------------------------------------------------
  // constructor & destructor
  // ---------------------------------------------------------------------------
 public:
  explicit TlDenseVector_ImplEigen(
      const TlDenseVectorObject::index_type size = 0);
  TlDenseVector_ImplEigen(const TlDenseVector_ImplEigen& rhs);
  TlDenseVector_ImplEigen(const std::vector<double>& rhs);
  TlDenseVector_ImplEigen(const double* p,
                          const TlDenseVectorObject::index_type size);
  TlDenseVector_ImplEigen(const VectorDataType& rhs);

  operator std::vector<double>() const;

  virtual ~TlDenseVector_ImplEigen();

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
  TlDenseVector_ImplEigen& operator=(const TlDenseVector_ImplEigen& rhs);

  TlDenseVector_ImplEigen& operator+=(const TlDenseVector_ImplEigen& rhs);
  TlDenseVector_ImplEigen& operator-=(const TlDenseVector_ImplEigen& rhs);
  TlDenseVector_ImplEigen& operator*=(const double rhs);
  TlDenseVector_ImplEigen& operator/=(const double rhs);

  // ---------------------------------------------------------------------------
  // operations
  // ---------------------------------------------------------------------------
 public:
  virtual double sum() const;
  virtual void sortByGreater();
  TlDenseVector_ImplEigen& dotInPlace(const TlDenseVector_ImplEigen& rhs);

  // ---------------------------------------------------------------------------
  // variables
  // ---------------------------------------------------------------------------
 protected:
  VectorDataType vector_;

  // ---------------------------------------------------------------------------
  // others
  // ---------------------------------------------------------------------------
  friend class TlDenseVector_ImplViennaCL;

  friend TlDenseVector_ImplEigen operator+(const TlDenseVector_ImplEigen& rhs1,
                                           const TlDenseVector_ImplEigen& rhs2);
  friend TlDenseVector_ImplEigen operator-(const TlDenseVector_ImplEigen& rhs1,
                                           const TlDenseVector_ImplEigen& rhs2);

  friend TlDenseVector_ImplEigen operator*(const TlDenseVector_ImplEigen& rhs1,
                                           const double rhs2);
  friend TlDenseVector_ImplEigen operator*(const double rhs1,
                                           const TlDenseVector_ImplEigen& rhs2);

  // DM(G) * DV
  friend TlDenseVector_ImplEigen operator*(
      const TlDenseGeneralMatrix_ImplEigen& mat,
      const TlDenseVector_ImplEigen& vec);
  // DV * DM(G)
  friend TlDenseVector_ImplEigen operator*(
      const TlDenseVector_ImplEigen& vec,
      const TlDenseGeneralMatrix_ImplEigen& mat);
  // DM(S) * DV
  friend TlDenseVector_ImplEigen operator*(
      const TlDenseSymmetricMatrix_ImplEigen& mat,
      const TlDenseVector_ImplEigen& vec);
  // DV * DM(S)
  friend TlDenseVector_ImplEigen operator*(
      const TlDenseVector_ImplEigen& vec,
      const TlDenseSymmetricMatrix_ImplEigen& mat);

  // SM(G) x DV
  friend TlDenseVector_ImplEigen operator*(const TlSparseGeneralMatrix_ImplEigen& mat, const TlDenseVector_ImplEigen& vtr);
  // DV x SM(G)
  friend TlDenseVector_ImplEigen operator*(const TlDenseVector_ImplEigen& vtr, const TlSparseGeneralMatrix_ImplEigen& mat);
  // SM(S) x DV
  friend TlDenseVector_ImplEigen operator*(const TlSparseSymmetricMatrix_ImplEigen& mat, const TlDenseVector_ImplEigen& vtr);

 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW  // Eigen macro
};

#endif  // TL_DENSE_VECTOR_IMPL_EIGEN_H
