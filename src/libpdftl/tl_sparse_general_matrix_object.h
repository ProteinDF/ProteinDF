#ifndef TL_SPARSE_GENERAL_MATRIX_OBJECT_H
#define TL_SPARSE_GENERAL_MATRIX_OBJECT_H

#include "tl_sparse_matrix_impl_object.h"
#include "tl_matrix_object.h"

class TlSparseGeneralMatrixObject : public TlMatrixObject {
    public:
  // ---------------------------------------------------------------------------
  // constructor & destructor
  // ---------------------------------------------------------------------------
  TlSparseGeneralMatrixObject(TlSparseMatrix_ImplObject* pImpl = NULL);
  virtual ~TlSparseGeneralMatrixObject();

 public:
  // ---------------------------------------------------------------------------
  // properties
  // ---------------------------------------------------------------------------
  virtual TlMatrixObject::index_type getNumOfRows() const;
  virtual TlMatrixObject::index_type getNumOfCols() const;
  virtual void resize(const TlMatrixObject::index_type row,
                      const TlMatrixObject::index_type col);

  virtual double get(const TlMatrixObject::index_type row,
                     const TlMatrixObject::index_type col) const;
  virtual void set(const TlMatrixObject::index_type row,
                   const TlMatrixObject::index_type col, const double value);
  virtual void add(const TlMatrixObject::index_type row,
                   const TlMatrixObject::index_type col, const double value);
  virtual void mul(const TlMatrixObject::index_type row,
                   const TlMatrixObject::index_type col, const double value);

  // ---------------------------------------------------------------------------
  // operators
  // ---------------------------------------------------------------------------

  // ---------------------------------------------------------------------------
  // I/O
  // ---------------------------------------------------------------------------
  virtual bool load(const std::string& filePath);
  virtual bool save(const std::string& filePath) const;

  // ---------------------------------------------------------------------------
  // variables
  // ---------------------------------------------------------------------------
 protected:
  TlSparseMatrix_ImplObject* pImpl_;
};

std::ostream& operator<<(std::ostream& stream,
                         const TlSparseGeneralMatrixObject& mat);

#endif // TL_SPARSE_GENERAL_MATRIX_OBJECT_H
