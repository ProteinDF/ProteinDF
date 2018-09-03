#ifndef TL_DENSE_SYMMETRIC_MATRIX_OBJECT_H
#define TL_DENSE_SYMMETRIC_MATRIX_OBJECT_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAVE_CONFIG_H

#include "TlSerializeData.h"
#include "tl_matrix_object.h"
class TlDenseMatrix_ImplObject;
class TlDenseGeneralMatrixObject;

class TlDenseSymmetricMatrixObject : public TlMatrixObject {
 public:
  // ---------------------------------------------------------------------------
  // constructor & destructor
  // ---------------------------------------------------------------------------
  TlDenseSymmetricMatrixObject(TlDenseMatrix_ImplObject* pImpl = NULL);
  virtual ~TlDenseSymmetricMatrixObject();

  virtual void vtr2mat(const std::vector<double>& vtr);
  
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
  virtual double getRMS() const;
  virtual double getMaxAbsoluteElement(
      TlMatrixObject::index_type* outRow = NULL,
      TlMatrixObject::index_type* outCol = NULL) const;

  void pivotedCholeskyDecomposition(TlDenseGeneralMatrixObject* pL,
                                    const double threshold) const;

  // ---------------------------------------------------------------------------
  // I/O
  // ---------------------------------------------------------------------------
  virtual bool load(const std::string& filePath);
  virtual bool save(const std::string& filePath) const;

  virtual bool saveText(const std::string& filePath) const;
  virtual void saveText(std::ostream& os) const;

#ifdef HAVE_HDF5
  virtual bool loadHdf5(const std::string& filepath, const std::string& h5path);
  virtual bool saveHdf5(const std::string& filepath,
                        const std::string& h5path) const;
#endif  // HAVE_HDF5

  void loadSerializeData(const TlSerializeData& data);
  TlSerializeData getSerializeData() const;

  // ---------------------------------------------------------------------------
  // variables
  // ---------------------------------------------------------------------------
 protected:
  TlDenseMatrix_ImplObject* pImpl_;
};

std::ostream& operator<<(std::ostream& stream,
                         const TlDenseSymmetricMatrixObject& mat);

#endif  // TL_DENSE_SYMMETRIC_MATRIX_OBJECT_H
