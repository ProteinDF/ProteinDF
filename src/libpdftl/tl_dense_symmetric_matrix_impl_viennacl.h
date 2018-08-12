#ifndef TL_DENSE_SYMMETRIC_MATRIX_IMPL_VIENNACL_H
#define TL_DENSE_SYMMETRIC_MATRIX_IMPL_VIENNACL_H

// IMPORTANT: Must be set prior to any ViennaCL includes if you want to use
// ViennaCL algorithms on Eigen objects
#define VIENNACL_WITH_EIGEN 1

#include "tl_dense_general_matrix_impl_viennacl.h"

class TlDenseGeneralMatrix_ImplViennaCL;
class TlDenseVector_ImplViennaCL;
class TlDenseSymmetricMatrix_ImplEigen;

class TlDenseSymmetricMatrix_ImplViennaCL
    : public TlDenseGeneralMatrix_ImplViennaCL {
  // ---------------------------------------------------------------------------
  // constructor & destructor
  // ---------------------------------------------------------------------------
 public:
  explicit TlDenseSymmetricMatrix_ImplViennaCL(
      const TlMatrixObject::index_type dim = 0);
  TlDenseSymmetricMatrix_ImplViennaCL(
      const TlDenseSymmetricMatrix_ImplViennaCL& rhs);
  TlDenseSymmetricMatrix_ImplViennaCL(
      const TlDenseGeneralMatrix_ImplViennaCL& rhs);
  TlDenseSymmetricMatrix_ImplViennaCL(
      const TlDenseSymmetricMatrix_ImplEigen& rhs);
  virtual ~TlDenseSymmetricMatrix_ImplViennaCL();

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
  // virtual double getMaxAbsoluteElement(TlMatrixObject::index_type* outRow,
  //                              TlMatrixObject::index_type* outCol) const;
  //
  // TlDenseSymmetetricMatrix_ImplViennaCL& dotInPlace(
  //     const TlDenseSymmetetricMatrix_ImplViennaCL& rhs);
  TlDenseSymmetricMatrix_ImplViennaCL transpose() const;
  TlDenseSymmetricMatrix_ImplViennaCL inverse() const;
  bool eig(TlDenseVector_ImplViennaCL* pEigval,
           TlDenseGeneralMatrix_ImplViennaCL* pEigvec) const;
  bool eig_QR(TlDenseVector_ImplViennaCL* pEigval,
              TlDenseGeneralMatrix_ImplViennaCL* pEigvec) const;

  // ---------------------------------------------------------------------------
  // private
  // ---------------------------------------------------------------------------

  // ---------------------------------------------------------------------------
  // others
  // ---------------------------------------------------------------------------
  friend class TlDenseGeneralMatrix_ImplViennaCL;
};

#endif  // TL_DENSE_SYMMETRIC_MATRIX_IMPL_VIENNACL_H
