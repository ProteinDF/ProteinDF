#ifndef TL_DENSE_GENERAL_MATRIX_IMPL_LAPACK_H
#define TL_DENSE_GENERAL_MATRIX_IMPL_LAPACK_H

#include "tl_dense_matrix_impl_object.h"

class TlDenseSymmetricMatrix_ImplLapack;
class TlDenseVector_ImplLapack;

class TlDenseGeneralMatrix_ImplLapack : public TlDenseMatrix_ImplObject {
 public:
  // ---------------------------------------------------------------------------
  // constructor & destructor
  // ---------------------------------------------------------------------------
 public:
  explicit TlDenseGeneralMatrix_ImplLapack(
      const TlMatrixObject::index_type row = 0,
      const TlMatrixObject::index_type col = 0);
  TlDenseGeneralMatrix_ImplLapack(const TlDenseGeneralMatrix_ImplLapack& rhs);
  TlDenseGeneralMatrix_ImplLapack(const TlDenseSymmetricMatrix_ImplLapack& rhs);
  virtual ~TlDenseGeneralMatrix_ImplLapack();

  // ---------------------------------------------------------------------------
  // static
  // ---------------------------------------------------------------------------
  static TlDenseGeneralMatrix_ImplLapack E(
      const TlMatrixObject::index_type dim);

  // ---------------------------------------------------------------------------
  // properties
  // ---------------------------------------------------------------------------
 public:
  virtual TlMatrixObject::index_type getNumOfRows() const;
  virtual TlMatrixObject::index_type getNumOfCols() const;
  virtual void resize(TlMatrixObject::index_type row,
                      TlMatrixObject::index_type col);

  virtual double get(const TlMatrixObject::index_type row,
                     const TlMatrixObject::index_type col) const;

  virtual void set(TlMatrixObject::index_type row,
                   TlMatrixObject::index_type col, double value);

  virtual void add(TlMatrixObject::index_type row,
                   TlMatrixObject::index_type col, double value);

  // ---------------------------------------------------------------------------
  // operators
  // ---------------------------------------------------------------------------
 public:
  TlDenseGeneralMatrix_ImplLapack& operator=(
      const TlDenseGeneralMatrix_ImplLapack& rhs);
  TlDenseGeneralMatrix_ImplLapack operator*(
      const TlDenseGeneralMatrix_ImplLapack& rhs) const;

  TlDenseGeneralMatrix_ImplLapack& operator+=(
      const TlDenseGeneralMatrix_ImplLapack& rhs);
  TlDenseGeneralMatrix_ImplLapack& operator-=(
      const TlDenseGeneralMatrix_ImplLapack& rhs);
  TlDenseGeneralMatrix_ImplLapack& operator*=(const double coef);
  TlDenseGeneralMatrix_ImplLapack& operator/=(const double coef);
  TlDenseGeneralMatrix_ImplLapack& operator*=(
      const TlDenseGeneralMatrix_ImplLapack& rhs);

  // ---------------------------------------------------------------------------
  // operations
  // ---------------------------------------------------------------------------
 public:
  // virtual double sum() const;
  // virtual double getRMS() const;
  // virtual double getMaxAbsoluteElement(
  //     TlMatrixObject::index_type* outRow,
  //     TlMatrixObject::index_type* outCol) const;

  const TlDenseGeneralMatrix_ImplLapack& dotInPlace(
      const TlDenseGeneralMatrix_ImplLapack& rhs);
  TlDenseGeneralMatrix_ImplLapack transpose() const;
  TlDenseGeneralMatrix_ImplLapack inverse() const;

  // ---------------------------------------------------------------------------
  // protected
  // ---------------------------------------------------------------------------
  virtual void initialize(bool clearIfNeeded = true);
  virtual TlMatrixObject::size_type getNumOfElements() const;
  virtual TlMatrixObject::size_type index(TlMatrixObject::index_type row,
                                          TlMatrixObject::index_type col) const;

  // ---------------------------------------------------------------------------
  // variables
  // ---------------------------------------------------------------------------
 protected:
  TlMatrixObject::index_type row_;
  TlMatrixObject::index_type col_;
  double* matrix_;  /// 行列要素

  // ---------------------------------------------------------------------------
  // friends
  // ---------------------------------------------------------------------------
  friend class TlDenseSymmetricMatrix_ImplLapack;
  friend TlDenseVector_ImplLapack operator*(
      const TlDenseGeneralMatrix_ImplLapack& mat,
      const TlDenseVector_ImplLapack& vec);
  friend TlDenseVector_ImplLapack operator*(
      const TlDenseVector_ImplLapack& vec,
      const TlDenseGeneralMatrix_ImplLapack& mat);
};

#endif  // TL_DENSE_GENERAL_MATRIX_IMPL_LAPACK_H
