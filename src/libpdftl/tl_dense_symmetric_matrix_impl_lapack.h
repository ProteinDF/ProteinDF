#ifndef TL_DENSE_SYMMETRIC_MATRIX_IMPL_LAPACK_H
#define TL_DENSE_SYMMETRIC_MATRIX_IMPL_LAPACK_H

#include "tl_dense_general_matrix_impl_lapack.h"
class TlDenseVector_ImplLapack;

class TlDenseSymmetricMatrix_ImplLapack
    : public TlDenseGeneralMatrix_ImplLapack {
  // ---------------------------------------------------------------------------
  // constructor & destructor
  // ---------------------------------------------------------------------------
 public:
  explicit TlDenseSymmetricMatrix_ImplLapack(
      const TlMatrixObject::index_type dim = 0);
  TlDenseSymmetricMatrix_ImplLapack(
      const TlDenseSymmetricMatrix_ImplLapack& rhs);
  TlDenseSymmetricMatrix_ImplLapack(const TlDenseGeneralMatrix_ImplLapack& rhs);
  virtual ~TlDenseSymmetricMatrix_ImplLapack();

  // ---------------------------------------------------------------------------
  // properties
  // ---------------------------------------------------------------------------
 public:
  virtual void resize(TlMatrixObject::index_type row,
                      TlMatrixObject::index_type col);

  // ---------------------------------------------------------------------------
  // operators
  // ---------------------------------------------------------------------------
 public:
  TlDenseSymmetricMatrix_ImplLapack& operator*=(const double coef);
  // const TlDenseGeneralMatrix_ImplLapack operator*=(
  //   const TlDenseSymmetricMatrix_ImplLapack& rhs);

  // ---------------------------------------------------------------------------
  // operations
  // ---------------------------------------------------------------------------
 public:
  // virtual double sum() const;
  // virtual double getRMS() const;
  // virtual double getMaxAbsoluteElement(
  //     TlMatrixObject::index_type* outRow,
  //     TlMatrixObject::index_type* outCol) const;

  // const TlDenseGeneralMatrix_ImplLapack& dotInPlace(
  //     const TlDenseGeneralMatrix_ImplLapack& rhs);
  TlDenseSymmetricMatrix_ImplLapack transpose() const;
  TlDenseSymmetricMatrix_ImplLapack inverse() const;

  bool eig(TlDenseVector_ImplLapack* pEigVal,
           TlDenseGeneralMatrix_ImplLapack* pEigVec) const;

  // ---------------------------------------------------------------------------
  // I/O
  // ---------------------------------------------------------------------------
 public:
  // virtual void dump(double* buf, const std::size_t size) const;
  // virtual void restore(const double* buf, const std::size_t size);

  // ---------------------------------------------------------------------------
  // protected
  // ---------------------------------------------------------------------------
 protected:
  virtual TlMatrixObject::size_type getNumOfElements() const;

  virtual TlMatrixObject::size_type index(TlMatrixObject::index_type row,
                                          TlMatrixObject::index_type col) const;

  // column-major
  // a11, a21, a31, a41, ..., am1, a12, a22, a32, ...                                          
  virtual void vtr2mat(const std::vector<double>& vtr);

  // ---------------------------------------------------------------------------
  // private
  // ---------------------------------------------------------------------------

  // ---------------------------------------------------------------------------
  // friends
  // ---------------------------------------------------------------------------
  friend class TlDenseSymmetricMatrix_Lapack;

  friend TlDenseGeneralMatrix_ImplLapack operator*(
      const TlDenseSymmetricMatrix_ImplLapack& rhs1,
      const TlDenseGeneralMatrix_ImplLapack& rhs2);
  friend TlDenseGeneralMatrix_ImplLapack operator*(
      const TlDenseGeneralMatrix_ImplLapack& rhs1,
      const TlDenseSymmetricMatrix_ImplLapack& rhs2);

  friend TlDenseVector_ImplLapack operator*(
      const TlDenseSymmetricMatrix_ImplLapack& mat,
      const TlDenseVector_ImplLapack& vec);
};

#endif  // TL_DENSE_SYMMETRIC_MATRIX_IMPL_LAPACK_H
