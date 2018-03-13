#ifndef TLMATRIX_MMAP_OBJECT_H
#define TLMATRIX_MMAP_OBJECT_H

#include <string>
#include "tl_matrix_object.h"

class TlDenseMatrixMmapObject : public TlMatrixObject {
 public:
  TlDenseMatrixMmapObject(const TlMatrixObject::MatrixType matrixType,
                          const std::string& filePath, const index_type row,
                          const index_type col);
  TlDenseMatrixMmapObject(const TlMatrixObject::MatrixType matrixType,
                          const std::string& filePath);
  virtual ~TlDenseMatrixMmapObject();

  void resize(const index_type newRow, const index_type newCol);

 public:
  virtual std::size_t getMemSize() const;

  virtual double get(index_type row, index_type col) const;
  virtual void set(index_type row, index_type col, double value);
  virtual void add(index_type row, index_type col, double value);

  virtual void setRowVector(const index_type row, const TlVector_BLAS& v);
  virtual void setColVector(const index_type col, const TlVector_BLAS& v);
  virtual TlVector_BLAS getRowVector(const index_type row) const;
  virtual TlVector_BLAS getColVector(const index_type col) const;

 protected:
  virtual TlDenseMatrixMmapObject* copy(const std::string& path) const = 0;
  virtual size_type getIndex(const index_type row,
                             const index_type col) const = 0;
  virtual TlMatrixObject::size_type getNumOfElements() const = 0;

 private:
  virtual bool load(const std::string& path);
  virtual bool save(const std::string& path) const;

 protected:
  void createNewFile();
  void openFile();

 private:
  void getHeaderInfo();
  void newMmap();
  void syncMmap();
  void deleteMmap();

 protected:
  std::string filePath_;

  char* mmapBegin_;
  double* dataBegin_;

  std::size_t headerSize_;
  std::size_t fileSize_;
};

#endif  // TLMATRIX_MMAP_OBJECT_H
