#ifndef TLMATRIX_RSFD_MMAP_H
#define TLMATRIX_RSFD_MMAP_H

#include "tl_dense_matrix_mmap_object.h"

class TlMatrixMmapRowMajor : public TlDenseMatrixMmapObject {
 public:
  TlMatrixMmapRowMajor(const std::string& filePath, const index_type row,
                       const index_type col);
  TlMatrixMmapRowMajor(const std::string& filePath);
  virtual ~TlMatrixMmapRowMajor();

 protected:
  virtual TlMatrixMmapRowMajor* copy(const std::string& path) const;
  virtual size_type getIndex(const index_type row, const index_type col) const;
  virtual TlMatrixObject::size_type getNumOfElements() const;
};

#endif  // TLMATRIX_RSFD_MMAP_H
