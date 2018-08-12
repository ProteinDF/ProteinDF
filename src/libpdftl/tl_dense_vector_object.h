#ifndef TL_DENSE_VECTOR_OBJECT_H
#define TL_DENSE_VECTOR_OBJECT_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAVE_CONFIG_H

#include <string>
#include "TlLogging.h"
#include "tl_vector_object.h"

class TlDenseVector_ImplObject;

class TlDenseVectorObject : public TlVectorObject {
  // ---------------------------------------------------------------------------
  // constructor & destructor
  // ---------------------------------------------------------------------------
 public:
  TlDenseVectorObject();
  virtual ~TlDenseVectorObject();

  // ---------------------------------------------------------------------------
  // properties
  // ---------------------------------------------------------------------------
 public:
  virtual size_type getSize() const;
  virtual void resize(const size_type newSize);

  virtual double get(const index_type index) const;
  virtual void set(const index_type index, const double value);
  virtual void add(const index_type index, const double value);
  virtual void mul(const index_type index, const double value);

  // ---------------------------------------------------------------------------
  // operators
  // ---------------------------------------------------------------------------

  // ---------------------------------------------------------------------------
  // operations
  // ---------------------------------------------------------------------------
  // virtual void push_back(const double value);
  virtual double getMaxAbsoluteElement(
      TlDenseVectorObject::index_type* index = NULL) const;
  virtual double sum() const;
  virtual double norm() const;
  virtual double norm2() const;
  virtual TlDenseVectorObject::index_type argmax(
      const TlDenseVectorObject::index_type& begin,
      const TlDenseVectorObject::index_type& end) const;

  // TlDenseVectorObject& dotInPlace(const TlVector_BLAS& rhs);
  void sortByGreater();

  // ---------------------------------------------------------------------------
  // I/O
  // ---------------------------------------------------------------------------
 public:
  virtual bool load(const std::string& filePath);
  virtual bool save(const std::string& filePath) const;
  virtual bool loadText(const std::string& filePath);

  virtual void outputText(std::ostream& os) const;

#ifdef HAVE_HDF5
  virtual bool saveHdf5(const std::string& filepath,
                        const std::string& h5path) const;
  virtual bool loadHdf5(const std::string& filepath, const std::string& h5path);
#endif  // HAVE_HDF5

  // ---------------------------------------------------------------------------
  // variables
  // ---------------------------------------------------------------------------
 protected:
  TlDenseVector_ImplObject* pImpl_;
  TlLogging& log_;
};

std::ostream& operator<<(std::ostream& stream, const TlDenseVectorObject& mat);

#endif  // TL_DENSE_VECTOR_OBJECT_H
