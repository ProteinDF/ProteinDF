#include <cassert>
#include <iomanip>
#include <iostream>
#include "TlFileMatrix.h"
#include "TlFile.h"

#define CACHE_GROUP_BIT (12)

TlFileMatrix::TlFileMatrix(const std::string& filePath, const int row, const int col, const size_t cacheSize)
        : filePath_(filePath), numOfRows_(row), numOfCols_(col), cacheCount_(0), cacheSize_(cacheSize)
{
    this->cache_.clear();
    if (this->cacheSize_ < 5 * 1024 * 1024) {
        this->cacheSize_ = 5 * 1024 * 1024;
    }
    this->open();
}

// called from sub-class only!
TlFileMatrix::TlFileMatrix(const std::string& filePath, const int row, const int col,
                           bool initialize, const size_t cacheSize)
        : filePath_(filePath), numOfRows_(row), numOfCols_(col), cacheCount_(0), cacheSize_(cacheSize)
{
    this->cache_.clear();
    if (this->cacheSize_ < 5 * 1024 * 1024) {
        this->cacheSize_ = 5 * 1024 * 1024;
    }
    if (initialize == true) {
        this->open();
    }
}

TlFileMatrix::~TlFileMatrix()
{
    std::list<CacheUnit>::const_iterator pEnd = this->cache_.end();
    for (std::list<CacheUnit>::const_iterator p = this->cache_.begin();
            p != pEnd; ++p) {
        if (p->isUpdate == true) {
            this->writeDisk(*p);
        }
    }

    this->fs_.flush();
    this->fs_.close();
}


std::size_t TlFileMatrix::getMemSize() const
{
    std::size_t answer = sizeof(TlFileMatrix) + this->cacheSize_;

    return answer;
}


void TlFileMatrix::open()
{
    if (TlFile::isExist(this->filePath_) == false) {
        // create new file
        this->fs_.open(this->filePath_.c_str(), std::ios::binary | std::ios::trunc | std::ios::out);

        // バッファ機能を停止する
        this->fs_ << std::setiosflags(std::ios::unitbuf);

        const int nType = 0; // means RSFD
        this->fs_.write(reinterpret_cast<const char*>(&nType), sizeof(int));
        this->fs_.write(reinterpret_cast<const char*>(&this->numOfRows_), sizeof(int));
        this->fs_.write(reinterpret_cast<const char*>(&this->numOfCols_), sizeof(int));

        const std::size_t size = std::size_t(this->numOfRows_) * std::size_t(this->numOfCols_);
        // size分ファイルにzeroを埋める
        const long unitSize = this->cacheSize_ / sizeof(double);
        ldiv_t ldivt = ldiv(size, unitSize);
        if (ldivt.quot > 0) {
            double* pBuf = new double[unitSize];
            for (int i = 0; i < unitSize; ++i) {
                pBuf[i] = 0.0;
            }

            for (long i = 0; i < ldivt.quot; ++i) {
                this->fs_.write(reinterpret_cast<const char*>(pBuf), sizeof(double) * unitSize);
            }
            delete[] pBuf;
            pBuf = NULL;
        }
        if (ldivt.rem > 0) {
            double* pBuf = new double[ldivt.rem];
            for (int i = 0; i < ldivt.rem; ++i) {
                pBuf[i] = 0.0;
            }

            this->fs_.write(reinterpret_cast<const char*>(pBuf), sizeof(double) * ldivt.rem);
            delete[] pBuf;
            pBuf = NULL;
        }

        // バッファ機能を再開する
        this->fs_ << std::resetiosflags(std::ios::unitbuf);
        this->fs_.flush();

        this->fs_.close();
    }

    this->fs_.open(this->filePath_.c_str(),
                   std::ios::binary | std::ios::in | std::ios::out);
    this->readHeader();
}

bool TlFileMatrix::readHeader()
{
    bool bAnswer = true;

    // get file size
    std::fstream::pos_type nFileSize = 0;
    {
        this->fs_.seekg(0, std::ios_base::beg);
        std::fstream::pos_type begin = this->fs_.tellg();
        this->fs_.seekg(0, std::ios_base::end);
        this->endPos_ = this->fs_.tellg();

        nFileSize = this->endPos_ - begin;
    }

    bool bCheckVariableType = false;
    this->startPos_ = 0;

    // int case:
    {
        int nType = 0;
        int nRows = 0;
        int nCols = 0;
        this->fs_.seekg(0, std::ios_base::beg);
        this->fs_.read((char*)&(nType), sizeof(int));
        this->fs_.read((char*)&(nRows), sizeof(int));
        this->fs_.read((char*)&(nCols), sizeof(int));

        const std::ifstream::pos_type nEstimatedFileSize =
            std::ifstream::pos_type(sizeof(int) *3)
            + (std::ifstream::pos_type(nRows) * std::ifstream::pos_type(nCols)
               * std::ifstream::pos_type(sizeof(double)));
        if (nEstimatedFileSize == nFileSize) {
            this->numOfRows_ = nRows;
            this->numOfCols_ = nCols;
            bCheckVariableType = true;
            this->startPos_ = sizeof(int) * 3;
        }
    }

    // long case:
    {
        int nType = 0;
        long nRows = 0;
        long nCols = 0;
        this->fs_.seekg(0, std::ios_base::beg);
        this->fs_.read((char*)&(nType), sizeof(int));
        this->fs_.read((char*)&(nRows), sizeof(long));
        this->fs_.read((char*)&(nCols), sizeof(long));

        const std::ifstream::pos_type nEstimatedFileSize =
            std::ifstream::pos_type(sizeof(int) + sizeof(long) * 2)
            + (std::ifstream::pos_type(nRows) * std::ifstream::pos_type(nCols)
               * std::ifstream::pos_type(sizeof(double)));
        if (nEstimatedFileSize == nFileSize) {
            this->numOfRows_ = nRows;
            this->numOfCols_ = nCols;
            bCheckVariableType = true;
            this->startPos_ = sizeof(int) + sizeof(long) * 2;
        }
    }

    if (bCheckVariableType == true) {
        this->fs_.seekg(this->startPos_, std::ios_base::beg);
    } else {
        std::cerr << "file size mismatch: " << this->filePath_ << std::endl;
        bAnswer = false;
    }

    return bAnswer;
}

void TlFileMatrix::set(const int row, const int col, const double value)
{
    *(this->getCachedData(row, col)) = value;
}

void TlFileMatrix::add(const int row, const int col, const double value)
{
    const double tmp = this->get(row, col) + value;
    *(this->getCachedData(row, col)) = tmp;
}

double TlFileMatrix::get(const int row, const int col) const
{
    return this->getCachedData(row, col);
}


TlVector TlFileMatrix::getRowVector(const int nRow) const
{
    assert((0 <= nRow) && (nRow < this->getNumOfCols()));

    const int nNumOfCols = this->getNumOfCols();
    TlVector answer(nNumOfCols);

    for (int i = 0; i < nNumOfCols; ++i) {
        answer[i] = this->get(nRow, i);
    }

    return answer;
}

TlVector TlFileMatrix::getColumnVector(const int nCol) const
{
    assert((0 <= nCol) && (nCol < this->getNumOfCols()));

    const int nNumOfRows = this->getNumOfRows();
    TlVector answer(nNumOfRows);

    for (int i = 0; i < nNumOfRows; ++i) {
        answer[i] = this->get(i, nCol);
    }

    return answer;
}

double* TlFileMatrix::getCachedData(const int row, const int col)
{
    static const unsigned long localIndexBit = (static_cast<unsigned long>(1) << CACHE_GROUP_BIT) -1;

    const size_t index = this->index(row, col);
    const size_t localIndex = (index & localIndexBit);
    this->updateCache(index);

    this->cache_.front().isUpdate = true;
    return &(this->cache_.front().data[localIndex]);
}

double TlFileMatrix::getCachedData(const int row, const int col) const
{
    static const unsigned long localIndexBit = (static_cast<unsigned long>(1) << CACHE_GROUP_BIT) -1;

    const size_t index = this->index(row, col);
    const size_t localIndex = (index & localIndexBit);
    this->updateCache(index);

    return this->cache_.front().data[localIndex];
}

void TlFileMatrix::updateCache(const size_t index) const
{
    // groupの探索
    const size_t group = (static_cast<unsigned long>(index) >> CACHE_GROUP_BIT);

    std::list<CacheUnit>::iterator pCacheUnitEnd = this->cache_.end();
    std::list<CacheUnit>::iterator found = std::find_if(this->cache_.begin(),
                                                        pCacheUnitEnd, CacheUnitComp(group));

    if (found != pCacheUnitEnd) {
        // 先頭に移動
        this->cache_.splice(this->cache_.begin(), this->cache_, found);
    } else {
        // groupの読み込み
        const size_t blockSize = sizeof(double) * (static_cast<unsigned long>(1) << CACHE_GROUP_BIT);
        const std::fstream::pos_type pos = std::fstream::pos_type(this->startPos_)
            + std::fstream::pos_type(group) * std::fstream::pos_type(blockSize);
        this->fs_.seekg(pos, std::ios_base::beg);

        const size_t groupCount = (1 << CACHE_GROUP_BIT);
        const size_t count = std::min(groupCount, (this->maxIndex() - group * groupCount));

        CacheUnit cu(group);
        cu.data.resize(count, 0.0);
        this->fs_.read((char*)&(cu.data[0]), sizeof(double) * count);

        this->cache_.push_front(cu);
        ++(this->cacheCount_);
    }

    // バッファサイズから溢れた分を削除
    static const size_t groupCacheSize = sizeof(double) * (static_cast<unsigned long>(1) << CACHE_GROUP_BIT);
    while ((this->cacheCount_ * groupCacheSize) > this->cacheSize_) {
        if (this->cache_.back().isUpdate == true) {
            this->writeDisk(this->cache_.back());
        }
        this->cache_.pop_back();
        --(this->cacheCount_);
    }
}

void TlFileMatrix::writeDisk(const TlFileMatrix::CacheUnit& cu) const
{
    std::string stateStr = "";
    std::ios_base::iostate state = this->fs_.rdstate();
    if ((state & std::ios_base::eofbit) || (state & std::ios_base::failbit)) {
        this->fs_.clear();
    }

    static const size_t groupCacheSize = sizeof(double) * (static_cast<unsigned long>(1) << CACHE_GROUP_BIT);
    const std::fstream::pos_type pos =
        this->startPos_ + std::fstream::pos_type(cu.group) * std::fstream::pos_type(groupCacheSize);
    this->fs_.seekp(pos, std::ios_base::beg);

    const size_t bufferSize = sizeof(double) * cu.data.size();

    this->fs_.write((char*)&(cu.data[0]), bufferSize);
}


TlFileMatrix& TlFileMatrix::operator*=(const double coef)
{
    // ToDo 高速化
    const index_type numOfRows = this->getNumOfRows();
    const index_type numOfCols = this->getNumOfCols();
    for (index_type r = 0; r < numOfRows; ++r) {
        for (index_type c = 0; c < numOfCols; ++c) {
            const double value = this->get(r, c) * coef;
            this->set(r, c, value);
        }
    }

    return *this;
}

TlMatrix TlFileMatrix::getBlockMatrix(const index_type row, const index_type col,
                                      const index_type rowDistance,
                                      const index_type colDistance) const
{
    assert((0 <= row) && (row < this->getNumOfRows()));
    assert((0 <= col) && (col < this->getNumOfCols()));
    assert(0 < rowDistance);
    assert(0 < colDistance);

    assert(0 <= (row + rowDistance) && (row + rowDistance) <= this->getNumOfRows());
    assert(0 <= (col + colDistance) && (col + colDistance) <= this->getNumOfCols());

    TlMatrix answer(rowDistance, colDistance);
#pragma omp parallel for
    for (index_type dr = 0; dr < rowDistance; ++dr) {
        const index_type r = row + dr;
        for (index_type dc = 0; dc < colDistance; ++dc) {
            const index_type c = col + dc;

            answer.set(dr, dc, this->get(r, c));
        }
    }

    return answer;
}

void TlFileMatrix::setBlockMatrix(const index_type row,
                                  const index_type col,
                                  const TlMatrix& matrix)
{
    const index_type row_distance = matrix.getNumOfRows();
    const index_type col_distance = matrix.getNumOfCols();

    assert(0 <= row && row < this->getNumOfRows());
    assert(0 <= col && col < this->getNumOfCols());
    assert(0 < (row + row_distance) && (row + row_distance) <= this->getNumOfRows());
    assert(0 < (col + col_distance) && (col + col_distance) <= this->getNumOfCols());

    // this->set()はスレッドセーフではないのでOpenMPでは注意すること
    for (index_type dr = 0; dr < row_distance; ++dr) {
        const index_type r = row + dr;
        for (index_type dc = 0; dc < col_distance; ++dc) {
            const index_type c = col + dc;

            this->set(r, c,  matrix.get(dr, dc));
        }
    }
}

