#ifndef TL_DENSE_GENERAL_MATRIX_OBJECT_H
#define TL_DENSE_GENERAL_MATRIX_OBJECT_H

#include "tl_dense_matrix_impl_object.h"
#include "tl_matrix_object.h"

class TlDenseGeneralMatrixObject : public TlMatrixObject {
 public:
  // ---------------------------------------------------------------------------
  // constructor & destructor
  // ---------------------------------------------------------------------------
  TlDenseGeneralMatrixObject(TlDenseMatrix_ImplObject* pImpl = NULL);
  virtual ~TlDenseGeneralMatrixObject();

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

  // ---------------------------------------------------------------------------
  // Operations
  // ---------------------------------------------------------------------------
  virtual double sum() const;
  virtual double trace() const;

  // ---------------------------------------------------------------------------
  // I/O
  // ---------------------------------------------------------------------------
  virtual bool load(const std::string& filePath);
  virtual bool save(const std::string& filePath) const;

#ifdef HAVE_HDF5
  virtual bool saveHdf5(const std::string& filepath,
                        const std::string& h5path) const;
  virtual bool loadHdf5(const std::string& filepath, const std::string& h5path);
#endif  // HAVE_HDF5

  // ---------------------------------------------------------------------------
  // variables
  // ---------------------------------------------------------------------------
 protected:
  TlDenseMatrix_ImplObject* pImpl_;
};

std::ostream& operator<<(std::ostream& stream,
                         const TlDenseGeneralMatrixObject& mat);

#endif  // TL_DENSE_GENERAL_MATRIX_OBJECT_H
