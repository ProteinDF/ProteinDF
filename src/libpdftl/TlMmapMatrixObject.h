#ifndef TLMMAPMATRIXOBJECT_H
#define TLMMAPMATRIXOBJECT_H

#include <string>
#include "TlMatrixObject.h"

class TlMmapMatrixObject : public TlMatrixObject {
 public:
  TlMmapMatrixObject(const TlMatrixObject::MatrixType matrixType,
                     const std::string& filePath, const index_type row,
                     const index_type col);
  virtual ~TlMmapMatrixObject();

  void resize(const index_type newRow, const index_type newCol);

 public:
  virtual index_type getNumOfRows() const;
  virtual index_type getNumOfCols() const;
  virtual std::size_t getMemSize() const;

  virtual double get(index_type row, index_type col) const;
  virtual void set(index_type row, index_type col, double value);
  virtual void add(index_type row, index_type col, double value);

  virtual void setRowVector(const index_type row, const TlVector& v);
  virtual void setColVector(const index_type col, const TlVector& v);
  virtual TlVector getRowVector(const index_type row) const;
  virtual TlVector getColVector(const index_type col) const;

 protected:
  virtual TlMmapMatrixObject* copy(const std::string& path) const = 0;
  virtual size_type getIndex(const index_type row,
                             const index_type col) const = 0;

 private:
  virtual bool load(const std::string& path);
  virtual bool save(const std::string& path) const;

 protected:
  void initialize();
  void createFile();
  void getHeaderInfo();
  void newMmap();
  void syncMmap();
  void deleteMmap();

 protected:
  index_type numOfRows_;
  index_type numOfCols_;

  std::string filePath_;

  char* mmapBegin_;
  double* dataBegin_;
  std::size_t fileSize_;
};

#endif  // TLMMAPMATRIXOBJECT_H
