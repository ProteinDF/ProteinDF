#include "tl_dense_general_matrix_mmap.h"

#include <cassert>
#include <iostream>

#include "TlFile.h"
#include "TlUtils.h"

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

// TlDenseGeneralMatrix_mmap* TlDenseGeneralMatrix_mmap::copy(
//     const std::string& path) const {
//     TlDenseGeneralMatrix_mmap* pObj = new TlDenseGeneralMatrix_mmap(path, 1,
//     1);

//     return pObj;
// }

TlMatrixObject::size_type TlDenseGeneralMatrix_mmap::getIndex(
    const index_type row, const index_type col) const {
    return TlMatrixObject::getIndex_CSFD(row, col);
}

TlMatrixObject::size_type TlDenseGeneralMatrix_mmap::getNumOfElements() const {
    return TlMatrixObject::getNumOfElements_CSFD();
}

void TlDenseGeneralMatrix_mmap::resize(const index_type newRow,
                                       const index_type newCol) {
    // backup
    const index_type oldRow = this->getNumOfRows();
    const index_type oldCol = this->getNumOfCols();
    this->deleteMmap();

    const std::string backupPath = this->filePath_ + ".bak";
    const int ret = TlFile::rename(this->filePath_, backupPath);
    assert(ret == 0);

    // new matrix
    this->row_ = newRow;
    this->col_ = newCol;
    this->createNewFile();
    this->openFile();
    assert(this->getNumOfRows() == newRow);
    assert(this->getNumOfCols() == newCol);

    TlDenseGeneralMatrix_mmap backup(backupPath);
    assert(oldRow == backup.getNumOfRows());
    assert(oldCol == backup.getNumOfCols());
    const index_type copyMaxRow = std::min(newRow, oldRow);
    const index_type copyMaxCol = std::min(newCol, oldCol);
    for (index_type c = 0; c < copyMaxCol; ++c) {
        for (index_type r = 0; r < copyMaxRow; ++r) {
            this->set(r, c, backup.get(r, c));
        }
    }

    TlFile::remove(backupPath);
}

void TlDenseGeneralMatrix_mmap::block(TlMatrixObject::index_type row,
                                      TlMatrixObject::index_type col,
                                      TlMatrixObject::index_type rowDistance,
                                      TlMatrixObject::index_type colDistance,
                                      TlDenseGeneralMatrixObject* pOut) const {
    assert((0 <= row) && (row < this->getNumOfRows()));
    assert((0 <= col) && (col < this->getNumOfCols()));
    assert(0 < rowDistance);
    assert(0 < colDistance);
    assert(0 <= (row + rowDistance) &&
           (row + rowDistance) <= this->getNumOfRows());
    assert(0 <= (col + colDistance) &&
           (col + colDistance) <= this->getNumOfCols());

    pOut->resize(rowDistance, colDistance);
#pragma omp parallel for
    for (TlMatrixObject::index_type dr = 0; dr < rowDistance; ++dr) {
        const TlMatrixObject::index_type r = row + dr;
        for (TlMatrixObject::index_type dc = 0; dc < colDistance; ++dc) {
            const TlMatrixObject::index_type c = col + dc;

            pOut->set(dr, dc, this->get(r, c));
        }
    }
}

void TlDenseGeneralMatrix_mmap::block(const TlMatrixObject::index_type row,
                                      const TlMatrixObject::index_type col,
                                      const TlDenseGeneralMatrixObject& ref) {
    const TlMatrixObject::index_type rowDistance = ref.getNumOfRows();
    const TlMatrixObject::index_type colDistance = ref.getNumOfCols();

    if (!((0 <= row && row < this->getNumOfRows()) &&
          (0 <= col && col < this->getNumOfCols()) &&
          (0 < (row + rowDistance) &&
           (row + rowDistance) <= this->getNumOfRows()) &&
          (0 < (col + colDistance) &&
           (col + colDistance) <= this->getNumOfCols()))) {
        this->log_.critical(TlUtils::format(
            "setBlockMatrix() start(%d, %d) mat(%d, %d) -> (%d, %d) @%s.%d",
            row, col, ref.getNumOfRows(), ref.getNumOfCols(),
            this->getNumOfRows(), this->getNumOfCols(), __FILE__, __LINE__));
    }
    assert(0 <= row && row < this->getNumOfRows());
    assert(0 <= col && col < this->getNumOfCols());
    assert(0 < (row + rowDistance) &&
           (row + rowDistance) <= this->getNumOfRows());
    assert(0 < (col + colDistance) &&
           (col + colDistance) <= this->getNumOfCols());

#pragma omp parallel for
    for (TlMatrixObject::index_type dr = 0; dr < rowDistance; ++dr) {
        const TlMatrixObject::index_type r = row + dr;
        for (TlMatrixObject::index_type dc = 0; dc < colDistance; ++dc) {
            const TlMatrixObject::index_type c = col + dc;

            this->set(r, c, ref.get(dr, dc));
        }
    }
}
