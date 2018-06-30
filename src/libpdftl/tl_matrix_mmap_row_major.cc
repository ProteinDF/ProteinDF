#include "tl_matrix_mmap_row_major.h"
#include <iostream>

TlMatrixMmapRowMajor::TlMatrixMmapRowMajor(const std::string& filePath,
                                           const index_type row,
                                           const index_type col)
    : TlDenseMatrixMmapObject(TlMatrixObject::RSFD, filePath, row, col) {
  this->createNewFile();
  this->openFile();
}

TlMatrixMmapRowMajor::TlMatrixMmapRowMajor(const std::string& filePath)
    : TlDenseMatrixMmapObject(TlMatrixObject::RSFD, filePath) {
  this->openFile();
}

TlMatrixMmapRowMajor::~TlMatrixMmapRowMajor() {}

TlMatrixMmapRowMajor* TlMatrixMmapRowMajor::copy(
    const std::string& path) const {
  TlMatrixMmapRowMajor* pObj = new TlMatrixMmapRowMajor(path, 1, 1);

  return pObj;
}

TlMatrixObject::size_type TlMatrixMmapRowMajor::getIndex(
    const index_type row, const index_type col) const {
  return TlMatrixObject::getIndex_RSFD(row, col);
}

TlMatrixObject::size_type TlMatrixMmapRowMajor::getNumOfElements() const {
  return TlMatrixObject::getNumOfElements_RSFD();
}
