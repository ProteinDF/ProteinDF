#ifndef TL_DENSE_MATRIX_IMPL_OBJECT_H
#define TL_DENSE_MATRIX_IMPL_OBJECT_H

#include <vector>

#include "TlLogging.h"
#include "tl_matrix_object.h"

class TlDenseVectorObject;

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

  /// Requests that the matrix capacity be at least enough to contain "row x col" elements.
  virtual void reserve(const TlMatrixObject::index_type row,
                       const TlMatrixObject::index_type col);

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
  virtual std::vector<double> diagonals() const;
  virtual double sum() const;
  virtual double trace() const;
  virtual double getRMS() const;
  virtual double getMaxAbsoluteElement(
      TlMatrixObject::index_type* outRow,
      TlMatrixObject::index_type* outCol) const;
  virtual void transposeInPlace() = 0;
  // const TlDenseGeneralMatrix_ImplEigen& dotInPlace(
  //     const TlDenseGeneralMatrix_ImplEigen& rhs);
  // TlDenseGeneralMatrix_ImplEigen inverse() const;

  // ---------------------------------------------------------------------------
  // I/O
  // ---------------------------------------------------------------------------
 public:
  void dump(double* buf, const std::size_t size) const;
  void restore(const double* buf, const std::size_t size);

  // ---------------------------------------------------------------------------
  // variables
  // ---------------------------------------------------------------------------
 protected:
  TlLogging& log_;
};

std::ostream& operator<<(std::ostream& stream,
                         const TlDenseMatrix_ImplObject& mat);

#endif  // TL_DENSE_MATRIX_IMPL_OBJECT_H
