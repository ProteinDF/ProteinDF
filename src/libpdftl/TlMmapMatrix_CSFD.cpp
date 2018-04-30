#include "TlMmapMatrix_CSFD.h"

TlMmapMatrix_CSFD::TlMmapMatrix_CSFD(const std::string& filePath,
                                     const index_type row, const index_type col)
    : TlMmapMatrixObject(TlMatrixObject::CSFD, filePath, row, col) {}

TlMmapMatrix_CSFD::~TlMmapMatrix_CSFD() {}

TlMmapMatrix_CSFD* TlMmapMatrix_CSFD::copy(const std::string& path) const {
  TlMmapMatrix_CSFD* pObj = new TlMmapMatrix_CSFD(path, 1, 1);

  return pObj;
}

TlMatrixObject::size_type TlMmapMatrix_CSFD::getIndex(
    const index_type row, const index_type col) const {
  return TlMatrixObject::getIndex_CSFD(row, col);
}
