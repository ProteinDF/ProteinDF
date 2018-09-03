#ifndef TL_SPARSE_MATRIX_IMPL_OBJECT_H
#define TL_SPARSE_MATRIX_IMPL_OBJECT_H

#include <vector>

#include "TlLogging.h"
#include "tl_matrix_object.h"

class TlSparseMatrix_ImplObject {
  // ---------------------------------------------------------------------------
  // constructor & destructor
  // ---------------------------------------------------------------------------
 public:
  TlSparseMatrix_ImplObject() : log_(TlLogging::getInstance()){};
  virtual ~TlSparseMatrix_ImplObject(){};

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
  virtual void mul(const TlMatrixObject::index_type row,
                   const TlMatrixObject::index_type col,
                   const double value) = 0;

  // ---------------------------------------------------------------------------
  // operations
  // ---------------------------------------------------------------------------
 public:

  // ---------------------------------------------------------------------------
  // I/O
  // ---------------------------------------------------------------------------
 public:

  // ---------------------------------------------------------------------------
  // variables
  // ---------------------------------------------------------------------------
 protected:
  TlLogging& log_;
};

std::ostream& operator<<(std::ostream& stream,
                         const TlSparseMatrix_ImplObject& mat);

#endif // TL_SPARSE_MATRIX_IMPL_OBJECT_H
