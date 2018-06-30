#ifndef TL_VECTOR_EIGEN_H
#define TL_VECTOR_EIGEN_H

#include <Eigen/Core>
#include <vector>
#include "tl_vector_abstract.h"

class TlEigenVector : public TlVectorAbstract {
 public:
  typedef Eigen::VectorXd VectorType;  // Matrix<double,Dynamic,1>
  typedef Eigen::Map<VectorType> MapType;
  typedef Eigen::Map<const VectorType> MapTypeConst;

 public:
  explicit TlEigenVector(TlVectorAbstract::index_type size = 0);
  TlEigenVector(const double* p, const TlVectorAbstract::size_type length);
  TlEigenVector(const TlEigenVector& rhs);
  TlEigenVector(const std::vector<double>& rhs);
  TlEigenVector(const VectorType& rhs);

  virtual ~TlEigenVector();

 public:
  TlVectorAbstract::size_type getSize() const;
  void resize(TlVectorAbstract::index_type newSize);

 public:
  double get(const TlVectorAbstract::index_type i) const;
  void set(const TlVectorAbstract::index_type i, const double value);
  void add(const TlVectorAbstract::index_type i, const double value);

 public:
  double operator[](const TlVectorAbstract::index_type i) const;
  double& operator[](const TlVectorAbstract::index_type i);

 public:
  bool load(const std::string& filePath);
  bool save(const std::string& filePath) const;

 protected:
  VectorType vector_;

 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW  // Eigen macro
};

#endif  // TL_VECTOR_EIGEN_H
