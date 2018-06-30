#ifndef TL_DENSE_SYMMETRIC_MATRIX_OBJECT_H
#define TL_DENSE_SYMMETRIC_MATRIX_OBJECT_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include "tl_dense_matrix_impl_object.h"
#include "tl_matrix_object.h"

class TlDenseSymmetricMatrixObject : public TlMatrixObject {
 public:
  // ---------------------------------------------------------------------------
  // constructor & destructor
  // ---------------------------------------------------------------------------
  TlDenseSymmetricMatrixObject(TlDenseMatrix_ImplObject* pImpl = NULL);
  virtual ~TlDenseSymmetricMatrixObject();

 public:
  // ---------------------------------------------------------------------------
  // properties
  // ---------------------------------------------------------------------------
  virtual TlMatrixObject::index_type getNumOfRows() const;
  virtual TlMatrixObject::index_type getNumOfCols() const;
  virtual void resize(const TlMatrixObject::index_type dim);

  virtual double get(const TlMatrixObject::index_type row,
                     const TlMatrixObject::index_type col) const;
  virtual void set(const TlMatrixObject::index_type row,
                   const TlMatrixObject::index_type col, const double value);
  virtual void add(const TlMatrixObject::index_type row,
                   const TlMatrixObject::index_type col, const double value);

  // ---------------------------------------------------------------------------
  // operations
  // ---------------------------------------------------------------------------
  virtual double sum() const;
  virtual double trace() const;

  // ---------------------------------------------------------------------------
  // I/O
  // ---------------------------------------------------------------------------
  virtual bool load(const std::string& filePath);
  virtual bool save(const std::string& filePath) const;

#ifdef HAVE_HDF5
  virtual bool loadHdf5(const std::string& filepath, const std::string& h5path);
  virtual bool saveHdf5(const std::string& filepath,
                        const std::string& h5path) const;
#endif  // HAVE_HDF5

  // ---------------------------------------------------------------------------
  // variables
  // ---------------------------------------------------------------------------
 protected:
  TlDenseMatrix_ImplObject* pImpl_;
};

std::ostream& operator<<(std::ostream& stream,
                         const TlDenseSymmetricMatrixObject& mat);

#endif  // TL_DENSE_SYMMETRIC_MATRIX_OBJECT_H
