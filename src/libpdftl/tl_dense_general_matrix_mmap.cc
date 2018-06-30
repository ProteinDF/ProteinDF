#include "tl_dense_general_matrix_mmap.h"

TlDenseGeneralMatrix_mmap::TlDenseGeneralMatrix_mmap(
    const std::string& filePath, const index_type row, const index_type col)
    : TlDenseMatrixMmapObject(TlMatrixObject::CSFD, filePath, row, col) {
  this->createNewFile();
  this->openFile();
}

TlDenseGeneralMatrix_mmap::TlDenseGeneralMatrix_mmap(
    const std::string& filePath)
    : TlDenseMatrixMmapObject(TlMatrixObject::CSFD, filePath) {
  this->openFile();
}

TlDenseGeneralMatrix_mmap::~TlDenseGeneralMatrix_mmap() {}

TlDenseGeneralMatrix_mmap* TlDenseGeneralMatrix_mmap::copy(
    const std::string& path) const {
  TlDenseGeneralMatrix_mmap* pObj = new TlDenseGeneralMatrix_mmap(path, 1, 1);

  return pObj;
}

TlMatrixObject::size_type TlDenseGeneralMatrix_mmap::getIndex(
    const index_type row, const index_type col) const {
  return TlMatrixObject::getIndex_CSFD(row, col);
}

TlMatrixObject::size_type TlDenseGeneralMatrix_mmap::getNumOfElements() const {
  return TlMatrixObject::getNumOfElements_CSFD();
}
