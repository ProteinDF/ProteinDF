#include "TlMmapMatrix_RSFD.h"

TlMmapMatrix_RSFD::TlMmapMatrix_RSFD(const std::string& filePath,
                                     const index_type row, const index_type col)
    : TlMmapMatrixObject(TlMatrixObject::RSFD, filePath, row, col) {}

TlMmapMatrix_RSFD::~TlMmapMatrix_RSFD() {}

TlMmapMatrix_RSFD* TlMmapMatrix_RSFD::copy(const std::string& path) const {
  TlMmapMatrix_RSFD* pObj = new TlMmapMatrix_RSFD(path, 1, 1);

  return pObj;
}

TlMatrixObject::size_type TlMmapMatrix_RSFD::getIndex(
    const index_type row, const index_type col) const {
  return TlMatrixObject::getIndex_RSFD(row, col);
}
