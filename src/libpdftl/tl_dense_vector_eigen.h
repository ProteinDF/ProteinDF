#ifndef TL_DENSE_VECTOR_EIGEN_H
#define TL_DENSE_VECTOR_EIGEN_H

#include <vector>

#include "tl_dense_vector_object.h"
class TlDenseGeneralMatrix_Eigen;
class TlDenseSymmetricMatrix_Eigen;
class TlDenseVector_ImplEigen;
class TlDenseVector_ViennaCL;
class TlSparseGeneralMatrix_Eigen;
class TlSparseSymmetricMatrix_Eigen;

class TlDenseVector_Eigen : public TlDenseVectorObject {
 public:
  explicit TlDenseVector_Eigen(TlDenseVectorObject::index_type size = 0);
  TlDenseVector_Eigen(const TlDenseVector_Eigen& rhs);
  TlDenseVector_Eigen(const std::vector<double>& rhs);
//   TlDenseVector_Eigen(const double* p,
//                       const TlDenseVectorObject::size_type length);
  TlDenseVector_Eigen(const TlDenseVector_ImplEigen& rhs);

#ifdef HAVE_VIENNACL
  TlDenseVector_Eigen(const TlDenseVector_ViennaCL& rhs);
#endif // HAVE_VIENNACL

  operator std::vector<double>() const;

  virtual ~TlDenseVector_Eigen();

  // ---------------------------------------------------------------------------
  // operators
  // ---------------------------------------------------------------------------
 public:
  TlDenseVector_Eigen& operator=(const TlDenseVector_Eigen& rhs);

  TlDenseVector_Eigen& operator+=(const TlDenseVector_Eigen& rhs);
  TlDenseVector_Eigen& operator-=(const TlDenseVector_Eigen& rhs);
  TlDenseVector_Eigen& operator*=(const double rhs);
  TlDenseVector_Eigen& operator/=(const double rhs);

  double operator*(const TlDenseVector_Eigen& rhs) const;

  // ---------------------------------------------------------------------------
  // operations
  // ---------------------------------------------------------------------------
  double dot(const TlDenseVector_Eigen& rhs) const;
  TlDenseVector_Eigen& dotInPlace(const TlDenseVector_Eigen& rhs);

  // ---------------------------------------------------------------------------
  // others
  // ---------------------------------------------------------------------------
  friend class TlDenseGeneralMatrix_Eigen;
  friend class TlDenseSymmetricMatrix_Eigen;
  friend class TlDenseVector_ViennaCL;

  friend TlDenseVector_Eigen operator+(const TlDenseVector_Eigen& rhs1,
                                       const TlDenseVector_Eigen& rhs2);
  friend TlDenseVector_Eigen operator-(const TlDenseVector_Eigen& rhs1,
                                       const TlDenseVector_Eigen& rhs2);
  friend TlDenseVector_Eigen operator*(const TlDenseVector_Eigen& rhs1,
                                       const double rhs2);
  friend TlDenseVector_Eigen operator*(const double rhs1,
                                       const TlDenseVector_Eigen& rhs2);

  // DM(G) * DV
  friend TlDenseVector_Eigen operator*(
      const TlDenseGeneralMatrix_Eigen& rhs1,
      const TlDenseVector_Eigen& rhs2);
  // DV * DM(G)
  friend TlDenseVector_Eigen operator*(
      const TlDenseVector_Eigen& rhs1,
      const TlDenseGeneralMatrix_Eigen& rhs2);
  // DM(S) * DV
  friend TlDenseVector_Eigen operator*(
      const TlDenseSymmetricMatrix_Eigen& rhs1,
      const TlDenseVector_Eigen& rhs2);
  // DV * DM(S)
  friend TlDenseVector_Eigen operator*(
      const TlDenseVector_Eigen& rhs1,
      const TlDenseSymmetricMatrix_Eigen& rhs2);

  friend TlDenseVector_Eigen operator*(const TlSparseGeneralMatrix_Eigen& mat, const TlDenseVector_Eigen& vtr);
  friend TlDenseVector_Eigen operator*(const TlSparseSymmetricMatrix_Eigen& mat, const TlDenseVector_Eigen& vtr);
};

#endif  // TL_DENSE_VECTOR_EIGEN_H
