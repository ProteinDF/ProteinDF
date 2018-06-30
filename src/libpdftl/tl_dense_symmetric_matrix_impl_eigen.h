#ifndef TL_DENSE_SYMMETRIC_MATRIX_IMPL_EIGEN_H
#define TL_DENSE_SYMMETRIC_MATRIX_IMPL_EIGEN_H

#include "tl_dense_general_matrix_impl_eigen.h"

class TlDenseVector_ImplEigen;

class TlDenseSymmetricMatrix_ImplEigen : public TlDenseGeneralMatrix_ImplEigen {
  // ---------------------------------------------------------------------------
  // constructor & destructor
  // ---------------------------------------------------------------------------
 public:
  explicit TlDenseSymmetricMatrix_ImplEigen(
      const TlMatrixObject::index_type dim = 0);
  TlDenseSymmetricMatrix_ImplEigen(const TlDenseSymmetricMatrix_ImplEigen& rhs);
  TlDenseSymmetricMatrix_ImplEigen(const TlDenseGeneralMatrix_ImplEigen& rhs);
  virtual ~TlDenseSymmetricMatrix_ImplEigen();

  // ---------------------------------------------------------------------------
  // properties
  // ---------------------------------------------------------------------------
  virtual void set(const TlMatrixObject::index_type row,
                   const TlMatrixObject::index_type col, const double value);

  virtual void add(const TlMatrixObject::index_type row,
                   const TlMatrixObject::index_type col, const double value);

  // ---------------------------------------------------------------------------
  // operators
  // ---------------------------------------------------------------------------

  // ---------------------------------------------------------------------------
  // operations
  // ---------------------------------------------------------------------------
 public:
  // virtual double sum() const;
  // virtual double getRMS() const;
  // virtual double getMaxAbsoluteElement(
  //     TlMatrixObject::index_type* outRow,
  //     TlMatrixObject::index_type* outCol) const;
  //
  // const TlDenseSymmetricMatrixMatrix_Eigen_new& dotInPlace(
  //     const TlDenseSymmetricMatrix_Eigen& rhs);
  TlDenseSymmetricMatrix_ImplEigen transpose() const;
  TlDenseSymmetricMatrix_ImplEigen inverse() const;
  bool eig(TlDenseVector_ImplEigen* pEigVal,
           TlDenseGeneralMatrix_ImplEigen* pEigVec) const;

  // ---------------------------------------------------------------------------
  // private
  // ---------------------------------------------------------------------------

  // ---------------------------------------------------------------------------
  // others
  // ---------------------------------------------------------------------------
  friend class TlDenseGeneralMatrix_ImplEigen;

 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW  // Eigen macro
};

#endif  // TL_DENSE_SYMMETRIC_MATRIX_IMPL_EIGEN_H
