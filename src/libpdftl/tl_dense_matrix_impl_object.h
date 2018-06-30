#ifndef TL_DENSE_MATRIX_IMPL_OBJECT_H
#define TL_DENSE_MATRIX_IMPL_OBJECT_H

#include "TlLogging.h"
#include "tl_matrix_object.h"

class TlDenseMatrix_ImplObject {
  // ---------------------------------------------------------------------------
  // constructor & destructor
  // ---------------------------------------------------------------------------
 public:
  TlDenseMatrix_ImplObject() : log_(TlLogging::getInstance()){};
  virtual ~TlDenseMatrix_ImplObject(){};

  // ---------------------------------------------------------------------------
  // properties
  // ---------------------------------------------------------------------------
 public:
  virtual TlMatrixObject::index_type getNumOfRows() const = 0;
  virtual TlMatrixObject::index_type getNumOfCols() const = 0;
  virtual void resize(const TlMatrixObject::index_type row,
                      const TlMatrixObject::index_type col) = 0;

  virtual double get(const TlMatrixObject::index_type row,
                     const TlMatrixObject::index_type col) const = 0;
  virtual void set(const TlMatrixObject::index_type row,
                   const TlMatrixObject::index_type col,
                   const double value) = 0;
  virtual void add(const TlMatrixObject::index_type row,
                   const TlMatrixObject::index_type col,
                   const double value) = 0;

  // ---------------------------------------------------------------------------
  // operations
  // ---------------------------------------------------------------------------
 public:
  virtual double sum() const;
  virtual double trace() const;
  virtual double getRMS() const;
  virtual double getMaxAbsoluteElement(
      TlMatrixObject::index_type* outRow,
      TlMatrixObject::index_type* outCol) const;

  // const TlDenseGeneralMatrix_ImplEigen& dotInPlace(
  //     const TlDenseGeneralMatrix_ImplEigen& rhs);
  // TlDenseGeneralMatrix_ImplEigen transpose() const;
  // TlDenseGeneralMatrix_ImplEigen inverse() const;

  // ---------------------------------------------------------------------------
  // variables
  // ---------------------------------------------------------------------------
 protected:
  TlLogging& log_;
};

#endif  // TL_DENSE_MATRIX_IMPL_OBJECT_H
