#include <errno.h>
#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "TlFile.h"
#include "TlUtils.h"
#include "tl_dense_matrix_mmap_object.h"
#include "tl_matrix_utils.h"

#define CLEAR_BUFSIZE (4096)

TlDenseMatrixMmapObject::TlDenseMatrixMmapObject(
    const TlMatrixObject::MatrixType matrixType, const std::string& filePath,
    const index_type row, const index_type col)
    : TlMatrixObject(matrixType, row, col), filePath_(filePath) {}

TlDenseMatrixMmapObject::TlDenseMatrixMmapObject(
    const TlMatrixObject::MatrixType matrixType, const std::string& filePath)
    : TlMatrixObject(matrixType, 1, 1), filePath_(filePath) {}

TlDenseMatrixMmapObject::~TlDenseMatrixMmapObject() {
    this->syncMmap();
    this->deleteMmap();
}

void TlDenseMatrixMmapObject::createNewFile() {
    if (TlFile::isExistFile(this->filePath_) == false) {
        // create new file
        std::fstream fs;
        fs.open(this->filePath_.c_str(),
                std::ios::binary | std::ios::trunc | std::ios::out);

        // バッファ機能を停止する
        fs << std::setiosflags(std::ios::unitbuf);

        const char matrixType = static_cast<char>(this->getType());
        fs.write(&matrixType, sizeof(char));
        fs.write(reinterpret_cast<const char*>(&this->row_),
                 sizeof(TlMatrixObject::index_type));
        fs.write(reinterpret_cast<const char*>(&this->col_),
                 sizeof(TlMatrixObject::index_type));

        const std::size_t size = this->getNumOfElements();
        // size分ファイルにzeroを埋める
        std::vector<double> buf(CLEAR_BUFSIZE);
        ldiv_t ldivt = ldiv(size, CLEAR_BUFSIZE);
        for (long i = 0; i < ldivt.quot; ++i) {
            fs.write(reinterpret_cast<const char*>(&(buf[0])),
                     sizeof(double) * CLEAR_BUFSIZE);
        }
        if (ldivt.rem > 0) {
            fs.write(reinterpret_cast<const char*>(&(buf[0])),
                     sizeof(double) * ldivt.rem);
        }

        // バッファ機能を再開する
        fs << std::resetiosflags(std::ios::unitbuf);
        fs.flush();

        fs.close();
    }
}

void TlDenseMatrixMmapObject::openFile() {
    this->getHeaderInfo();
    this->newMmap();
}

void TlDenseMatrixMmapObject::getHeaderInfo() {
    TlMatrixObject::MatrixType matrixType;
    index_type row = 0;
    index_type col = 0;

    TlMatrixUtils::FileSize headerSize =
        TlMatrixUtils::getHeaderInfo(this->filePath_, &matrixType, &row, &col);
    if ((headerSize > 0) && (matrixType == this->getType())) {
        this->row_ = row;
        this->col_ = col;
    } else {
        std::cerr << TlUtils::format("cannot read matrix header: %s@%d",
                                     __FILE__, __LINE__)
                  << std::endl;
        std::cerr << TlUtils::format("matrix: type=%d, row=%d, col=%d; hs=%ld",
                                     (int)matrixType, row, col, headerSize)
                  << std::endl;
    }

    this->headerSize_ = headerSize;
    this->fileSize_ = headerSize + sizeof(double) * this->getNumOfElements();
}

void TlDenseMatrixMmapObject::newMmap() {
    const int fd = open(this->filePath_.c_str(), O_RDWR, S_IRUSR | S_IWUSR);
    if (fd == -1) {
        int errorNo = errno;
        std::string errStr(strerror(errorNo));
        this->log_.critical(TlUtils::format(
            "open error: %s (%s)", this->filePath_.c_str(), errStr.c_str()));
        throw errStr;
    }

    struct stat sb;
    fstat(fd, &sb);
    assert(std::size_t(sb.st_size) == this->fileSize_);

    // const std::size_t pageSize = sysconf(_SC_PAGE_SIZE);
    // const std::size_t mapSize = (this->fileSize_ / pageSize +1) * pageSize;
    this->mmapBegin_ =
        (char*)mmap(NULL, this->fileSize_, PROT_READ | PROT_WRITE, MAP_SHARED,
                    fd, 0);  //  MAP_HUGETLB
    this->dataBegin_ = (double*)(this->mmapBegin_ + this->headerSize_);
    msync(this->mmapBegin_, this->getNumOfElements(), MS_ASYNC);
    madvise(this->mmapBegin_, this->getNumOfElements(), MADV_WILLNEED);

    // we can close the file after mmap() was called.
    close(fd);
}

void TlDenseMatrixMmapObject::syncMmap() {
    msync(this->mmapBegin_, this->fileSize_, MS_SYNC);
}

void TlDenseMatrixMmapObject::deleteMmap() {
    this->syncMmap();
    munmap(this->mmapBegin_, this->fileSize_);

    this->mmapBegin_ = NULL;
}

void TlDenseMatrixMmapObject::resize(const index_type newRow,
                                     const index_type newCol) {
    // backup
    const index_type oldRow = this->getNumOfRows();
    const index_type oldCol = this->getNumOfCols();
    this->syncMmap();
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

    TlDenseMatrixMmapObject* pBackup = this->copy(backupPath);
    assert(oldRow == pBackup->getNumOfRows());
    assert(oldCol == pBackup->getNumOfCols());
    const index_type copyMaxRow = std::min(newRow, oldRow);
    const index_type copyMaxCol = std::min(newCol, oldCol);
    for (index_type c = 0; c < copyMaxCol; ++c) {
        for (index_type r = 0; r < copyMaxRow; ++r) {
            this->set(r, c, pBackup->get(r, c));
        }
    }
    delete pBackup;
    pBackup = NULL;

    TlFile::remove(backupPath);
}

std::size_t TlDenseMatrixMmapObject::getMemSize() const { return 0; }

double TlDenseMatrixMmapObject::get(const index_type row,
                                    const index_type col) const {
    assert((0 <= row) && (row < this->getNumOfRows()));
    assert((0 <= col) && (col < this->getNumOfCols()));
    const size_type index = this->getIndex(row, col);

    return this->dataBegin_[index];
}

void TlDenseMatrixMmapObject::set(const index_type row, const index_type col,
                                  const double value) {
    assert((0 <= row) && (row < this->getNumOfRows()));
    assert((0 <= col) && (col < this->getNumOfCols()));

    const size_type index = this->getIndex(row, col);
    this->dataBegin_[index] = value;
}

void TlDenseMatrixMmapObject::add(const index_type row, const index_type col,
                                  const double value) {
    assert((0 <= row) && (row < this->getNumOfRows()));
    assert((0 <= col) && (col < this->getNumOfCols()));

    const size_type index = this->getIndex(row, col);
    this->dataBegin_[index] += value;
}

void TlDenseMatrixMmapObject::setRowVector(const index_type row,
                                           const TlDenseVector_Lapack& v) {
    const index_type numOfCols = this->getNumOfCols();
    assert(v.getSize() == numOfCols);

#pragma omp parallel for schedule(runtime)
    for (index_type i = 0; i < numOfCols; ++i) {
        this->set(row, i, v.get(i));
    }
}

void TlDenseMatrixMmapObject::setColVector(const index_type col,
                                           const TlDenseVector_Lapack& v) {
    const index_type numOfRows = this->getNumOfRows();
    assert(v.getSize() == numOfRows);

#pragma omp parallel for schedule(runtime)
    for (index_type i = 0; i < numOfRows; ++i) {
        this->set(i, col, v.get(i));
    }
}

std::vector<double> TlDenseMatrixMmapObject::getRowVector(
    const index_type row) const {
    assert((0 <= row) && (row < this->getNumOfRows()));

    const index_type numOfCols = this->getNumOfCols();
    std::vector<double> answer(numOfCols);

#pragma omp parallel for schedule(runtime)
    for (index_type i = 0; i < numOfCols; ++i) {
        answer[i] = this->get(row, i);
    }

    return answer;
}

std::vector<double> TlDenseMatrixMmapObject::getColVector(
    const index_type col) const {
    assert((0 <= col) && (col < this->getNumOfCols()));

    const index_type numOfRows = this->getNumOfRows();
    std::vector<double> answer(numOfRows);

#pragma omp parallel for schedule(runtime)
    for (index_type i = 0; i < numOfRows; ++i) {
        answer[i] = this->get(i, col);
    }

    return answer;
}

std::size_t TlDenseMatrixMmapObject::getRowVector(
    const index_type row, double* pBuf, const std::size_t count) const {
    assert((0 <= row) && (row < this->getNumOfRows()));

    const std::size_t numOfCols = this->getNumOfCols();
    const std::size_t copyCount = std::min(count, numOfCols);

#pragma omp parallel for schedule(runtime)
    for (std::size_t i = 0; i < copyCount; ++i) {
        pBuf[i] = this->get(row, i);
    }

    return copyCount;
}

std::size_t TlDenseMatrixMmapObject::getColVector(
    const index_type col, double* pBuf, const std::size_t count) const {
    assert((0 <= col) && (col < this->getNumOfCols()));

    const std::size_t numOfRows = this->getNumOfRows();
    const std::size_t copyCount = std::min(count, numOfRows);

#pragma omp parallel for schedule(runtime)
    for (index_type i = 0; i < copyCount; ++i) {
        pBuf[i] = this->get(i, col);
    }

    return copyCount;
}

bool TlDenseMatrixMmapObject::load(const std::string& path) { return true; }

bool TlDenseMatrixMmapObject::save(const std::string& path) const {
    return true;
}
