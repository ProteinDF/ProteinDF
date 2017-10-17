#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>

#include <iostream>
#include <fstream>
#include <iomanip>

#include "TlMmapMatrixObject.h"
#include "TlMatrix.h"
#include "TlFile.h"

#define CLEAR_BUFSIZE (4096)


TlMmapMatrixObject::TlMmapMatrixObject(const TlMatrixObject::MatrixType matrixType,
                                       const std::string& filePath,
                                       const index_type row, const index_type col)
    : TlMatrixObject(matrixType), numOfRows_(row), numOfCols_(col), filePath_(filePath)
{
    this->initialize();
}


TlMmapMatrixObject::~TlMmapMatrixObject()
{
    this->syncMmap();
    this->deleteMmap();
}


void TlMmapMatrixObject::initialize()
{
    this->createFile();
    this->getHeaderInfo();
    this->newMmap();
}


void TlMmapMatrixObject::createFile()
{
    if (TlFile::isExistFile(this->filePath_) == false) {
        // create new file
        std::fstream fs;
        fs.open(this->filePath_.c_str(), std::ios::binary | std::ios::trunc | std::ios::out);
        
        // バッファ機能を停止する
        fs << std::setiosflags(std::ios::unitbuf);
        
        const int matrixType = this->getType();
        fs.write(reinterpret_cast<const char*>(&matrixType), sizeof(int));
        fs.write(reinterpret_cast<const char*>(&this->numOfRows_), sizeof(int));
        fs.write(reinterpret_cast<const char*>(&this->numOfCols_), sizeof(int));
        
        const std::size_t size = std::size_t(this->numOfRows_) * std::size_t(this->numOfCols_);
        // size分ファイルにzeroを埋める
        std::vector<double> buf(CLEAR_BUFSIZE);
        ldiv_t ldivt = ldiv(size, CLEAR_BUFSIZE);
        for (long i = 0; i < ldivt.quot; ++i) {
            fs.write(reinterpret_cast<const char*>(&(buf[0])), sizeof(double) * CLEAR_BUFSIZE);
        }
        if (ldivt.rem > 0) {
            fs.write(reinterpret_cast<const char*>(&(buf[0])), sizeof(double) * ldivt.rem);
        }
        
        // バッファ機能を再開する
        fs << std::resetiosflags(std::ios::unitbuf);
        fs.flush();
        
        fs.close();
    }
}


void TlMmapMatrixObject::getHeaderInfo()
{
    int matrixType;
    index_type numOfRows = 0;
    index_type numOfCols = 0;
    if (TlMatrix::getHeaderInfo(this->filePath_, &matrixType, &numOfRows, &numOfCols)) {
        this->numOfRows_ = numOfRows;
        this->numOfCols_ = numOfCols;
    }
    
    this->fileSize_ = sizeof(int) * 3 + sizeof(double) * std::size_t(this->getNumOfRows()) * std::size_t(this->getNumOfCols());
}

        
void TlMmapMatrixObject::newMmap()
{
    const int fd = open(this->filePath_.c_str(), O_RDWR, S_IRUSR | S_IWUSR);
    if (fd == -1) {
        int errorNo = errno;
        std::string errStr(strerror(errorNo));
        this->log_.critical(TlUtils::format("open error: %s (%s)",
                                            this->filePath_.c_str(),
                                            errStr.c_str()));
        throw errStr;
    }

    struct stat sb;
    fstat(fd, &sb);
    assert(std::size_t(sb.st_size) == this->fileSize_);
    
    //const std::size_t pageSize = sysconf(_SC_PAGE_SIZE);
    //const std::size_t mapSize = (this->fileSize_ / pageSize +1) * pageSize;
    this->mmapBegin_ = (char*)mmap(NULL, this->fileSize_,  PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    this->dataBegin_ = (double*)(this->mmapBegin_ + sizeof(int) * 3);
    
    // we can close the file after mmap() was called.
    close(fd);
}


void TlMmapMatrixObject::syncMmap()
{
    msync(this->mmapBegin_, this->fileSize_, MS_SYNC);
}


void TlMmapMatrixObject::deleteMmap()
{
    this->syncMmap();
    munmap(this->mmapBegin_, this->fileSize_);

    this->mmapBegin_ = NULL;
}


void TlMmapMatrixObject::resize(const index_type newRow, const index_type newCol)
{
    // backup
    const index_type oldRow = this->getNumOfRows();
    const index_type oldCol = this->getNumOfCols();
    const std::string backupPath = TlFile::getTempFilePath();
    this->syncMmap();
    this->deleteMmap();
    TlFile::copy(this->filePath_, backupPath);

    // new matrix
    TlFile::remove(this->filePath_);
    this->numOfRows_ = newRow;
    this->numOfCols_ = newCol;
    this->initialize();
    assert(this->getNumOfRows() == newRow);
    assert(this->getNumOfCols() == newCol);
    
    // TlMmapMatrixObject backup(backupPath);
    TlMmapMatrixObject* pBackup = this->copy(backupPath);
    const index_type copyMaxRow = std::min(newRow, oldRow);
    const index_type copyMaxCol = std::min(newCol, oldCol);
    for (index_type c = 0; c < copyMaxCol; ++c) {
        for (index_type r = 0; r < copyMaxRow; ++r) {
            this->set(r, c, pBackup->get(r, c));
        }
    }
    delete pBackup;
    pBackup = NULL;

    // delete backup
    TlFile::remove(backupPath);
}


TlMatrixObject::index_type TlMmapMatrixObject::getNumOfRows() const
{
    return this->numOfRows_;
}


TlMatrixObject::index_type TlMmapMatrixObject::getNumOfCols() const
{
    return this->numOfCols_;
}


std::size_t TlMmapMatrixObject::getMemSize() const
{
    return 0;
}


double TlMmapMatrixObject::get(const index_type row, const index_type col) const
{
    assert((0 <= row) && (row < this->getNumOfRows()));
    assert((0 <= col) && (col < this->getNumOfCols()));
    const size_type index = this->getIndex_RSFD(row, col);

    return this->dataBegin_[index];
}


void TlMmapMatrixObject::set(const index_type row, const index_type col, const double value)
{
    assert((0 <= row) && (row < this->getNumOfRows()));
    assert((0 <= col) && (col < this->getNumOfCols()));
           
    const size_type index = this->getIndex_RSFD(row, col);
    this->dataBegin_[index] = value;
}


void TlMmapMatrixObject::add(const index_type row, const index_type col, const double value)
{
    assert((0 <= row) && (row < this->getNumOfRows()));
    assert((0 <= col) && (col < this->getNumOfCols()));
           
    const size_type index = this->getIndex_RSFD(row, col);
    this->dataBegin_[index] += value;
}


void TlMmapMatrixObject::setRowVector(const index_type row, const TlVector& v)
{
    const index_type numOfCols = this->getNumOfCols();
    assert(v.getSize() == numOfCols);
    
    for (index_type i = 0; i < numOfCols; ++i) {
        this->set(row, i, v[i]);
    }
}


void TlMmapMatrixObject::setColVector(const index_type col, const TlVector& v)
{
    const index_type numOfRows = this->getNumOfRows();
    assert(v.getSize() == numOfRows);

    for (index_type i = 0; i < numOfRows; ++i) {
        this->set(i, col, v[i]);
    }
}


TlVector TlMmapMatrixObject::getRowVector(const index_type row) const
{
    assert((0 <= row) && (row < this->getNumOfRows()));

    const index_type numOfCols = this->getNumOfCols();
    TlVector answer(numOfCols);

    for (index_type i = 0; i < numOfCols; ++i) {
        answer[i] = this->get(row, i);
    }

    return answer;
}


TlVector TlMmapMatrixObject::getColVector(const index_type col) const
{
    assert((0 <= col) && (col < this->getNumOfCols()));

    const index_type numOfRows = this->getNumOfRows();
    TlVector answer(numOfRows);

    for (index_type i = 0; i < numOfRows; ++i) {
        answer[i] = this->get(i, col);
    }

    return answer;
}


bool TlMmapMatrixObject::load(const std::string& path)
{
    return true;
}


bool TlMmapMatrixObject::save(const std::string& path) const
{
    return true;
}


