#ifdef HAVE_CONFIG_H
#include "config.h"    // this file created by autotools
#endif // HAVE_CONFIG_H

#include <iostream>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <numeric>
#include <limits>

#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"
#include "TlVector.h"
#include "TlUtils.h"
#include "TlMemManager.h"
#include "TlLogging.h"

const TlMatrix::size_type TlMatrix::MAX_LOOP = std::numeric_limits<int>::max();
bool TlMatrix::isUsingMemManagerDefault_ = false;

////////////////////////////////////////////////////////////////////////
//
void TlMatrix::useMemManager(bool isUsingMemManager)
{
    TlMatrix::isUsingMemManagerDefault_ = isUsingMemManager;
}

TlMatrix::TlMatrix(const index_type nRow, const index_type nCol)
    : TlMatrixObject(), m_nRows(nRow), m_nCols(nCol), data_(NULL),
      isUsingMemManager_(TlMatrix::isUsingMemManagerDefault_)
{
    assert((0 <= nRow) && (0 <= nCol));
    this->initialize();
}


// for sub-class
TlMatrix::TlMatrix(const index_type nRow, const index_type nCol, double* pData)
    : TlMatrixObject(), m_nRows(nRow), m_nCols(nCol), data_(pData),
      isUsingMemManager_(TlMatrix::isUsingMemManagerDefault_)
{
}


TlMatrix::TlMatrix(const TlMatrix& rhs)
    : m_nRows(rhs.getNumOfRows()), m_nCols(rhs.getNumOfCols()), data_(NULL),
      isUsingMemManager_(TlMatrix::isUsingMemManagerDefault_)
{

    this->initialize(false);
    const std::size_t size = this->getNumOfRows() * this->getNumOfCols();
    std::copy(rhs.data_, rhs.data_ + size, this->data_);
}


TlMatrix::TlMatrix(const TlSymmetricMatrix& rhs)
    : m_nRows(rhs.getNumOfRows()), m_nCols(rhs.getNumOfCols()), data_(NULL),
      isUsingMemManager_(TlMatrix::isUsingMemManagerDefault_)
{

    this->initialize(false);

    const int numOfRows = this->getNumOfRows();

    // FX1用富士通コンパイラでは
    // OpenMPをコンストラクタで使ってはいけない。
    for (int row = 0; row < numOfRows; ++row) {
        for (int col = 0; col < row ; ++col) {
            const double tmp = rhs.get(row, col);

            this->set(row, col, tmp);
            this->set(col, row, tmp);
        }

        // case : row == col
        this->set(row, row, rhs.get(row, row));
    }
}


TlMatrix::TlMatrix(const TlSerializeData& data)
    : m_nRows(1), m_nCols(1), data_(NULL),
      isUsingMemManager_(TlMatrix::isUsingMemManagerDefault_)
{
    this->m_nRows = std::max(data["row"].getInt(), 1);
    this->m_nCols = std::max(data["col"].getInt(), 1);
    this->initialize(false);

    const size_type size = this->getNumOfElements();
    for (size_type index = 0; index < size; ++index) {
        this->data_[index] = data["data"].getAt(index).getDouble();
    }
}


TlMatrix::TlMatrix(const TlVector& vct,
                   const index_type row, const index_type col)
    : TlMatrixObject(), m_nRows(row), m_nCols(col), data_(NULL),
      isUsingMemManager_(TlMatrix::isUsingMemManagerDefault_)
{
    assert((0 <= row) && (0 <= col));

    const size_type size = this->getNumOfElements();
    assert(vct.getSize() == size);

    this->initialize(false);

#ifndef __FUJITSU // cannot use OpenMP in constructor
#pragma omp parallel for
#endif // __FUJITSU 
    for (size_type index = 0; index < size; ++index) {
        this->data_[index] = vct[index];
    }
}


TlMatrix::~TlMatrix()
{
    this->clear();
}


int TlMatrix::getNumOfRows() const
{
    return this->m_nRows;
}


int TlMatrix::getNumOfCols() const
{
    return this->m_nCols;
}


void TlMatrix::initialize(bool isZeroClear)
{
#ifdef HAVE_MMAP
    if (this->isUsingMemManager_ == true) {
        this->initialize_usingMemManager(isZeroClear);
    } else {
        this->initialize_usingStandard(isZeroClear);
    }
#else
    this->initialize_usingStandard(isZeroClear);
#endif // HAVE_MMAP    
}


void TlMatrix::initialize_usingStandard(bool isZeroClear)
{
    const std::size_t size = this->getNumOfElements();
    this->data_ = new double[size];
    
    if (isZeroClear == true) {
        std::fill(this->data_, this->data_ + size, 0.0);
    }
}


void TlMatrix::initialize_usingMemManager(bool isZeroClear)
{
    TlMemManager& rMemManager = TlMemManager::getInstance();
    const std::size_t size = this->getNumOfElements();
    this->data_ = (double*)rMemManager.allocate(sizeof(double) * size);

    if (isZeroClear == true) {
        std::fill(this->data_, this->data_ + size, 0.0);
    }
}


void TlMatrix::clear()
{
    if (this->data_ != NULL) {
#ifdef HAVE_MMAP
        if (this->isUsingMemManager_ == true) {
            this->clear_usingMemManager();
        } else {
            this->clear_usingStandard();
        }
#else
        this->clear_usingStandard();
#endif // HAVE_MMAP
    }
}


void TlMatrix::clear_usingStandard()
{
    delete[] this->data_;
    this->data_ = NULL;
}


void TlMatrix::clear_usingMemManager()
{
    TlMemManager& rMemManager = TlMemManager::getInstance();
    rMemManager.deallocate((char*)this->data_, sizeof(double) * this->getNumOfElements());
    this->data_ = NULL;
}


void TlMatrix::resize(const int nRow, const int nCol)
{
    assert((nRow > 0) && (nCol > 0));

    TlMatrix oldMatrix(*this);

    this->clear();
    this->m_nRows = nRow;
    this->m_nCols = nCol;
    this->initialize(true);

    const int nMaxRowsForCopy = std::min(oldMatrix.getNumOfRows(), nRow);
    const int nMaxColsForCopy = std::min(oldMatrix.getNumOfCols(), nCol);
#pragma omp parallel for
    for (int c = 0; c < nMaxColsForCopy; ++c) {
        for (int r = 0; r < nMaxRowsForCopy; ++r) {
            this->set(r, c, oldMatrix.get(r, c));
        }
    }
}


std::size_t TlMatrix::index(index_type row,
                            index_type col) const
{
    assert((0 <= row) && (row < this->m_nRows));
    assert((0 <= col) && (col < this->m_nCols));

    return (row  + (this->m_nRows * col));
}


double TlMatrix::get(const index_type row,
                     const index_type col) const
{
    return this->data_[this->index(row, col)];
}


void TlMatrix::set(const index_type row,
                   const index_type col,
                   const double value)
{
    const std::size_t index = this->index(row, col);
    
#pragma omp critical(TlMatrix__set)
    {
        this->data_[index] = value;
    }
}


void TlMatrix::add(const index_type row,
                   const index_type col,
                   const double value)
{
    const std::size_t index = this->index(row, col);

#pragma omp atomic
    this->data_[index] += value;
}


TlVector TlMatrix::getVector() const
{
    return TlVector(this->data_, this->getNumOfElements());
}


TlMatrix TlMatrix::getBlockMatrix(const int nRow, const int nCol,
                                  const int nRowDistance, const int nColDistance) const
{
    assert((0 <= nRow) && (nRow < this->getNumOfRows()));
    assert((0 <= nCol) && (nCol < this->getNumOfCols()));
    assert(0 < nRowDistance);
    assert(0 < nColDistance);

    // for debug
    //   {
    //     std::cout << TlUtils::format("matrix size = (%d, %d)", this->getNumOfRows(), this->getNumOfCols())
    //          << std::endl;
    //     std::cout << TlUtils::format("(%d, %d) -> (%d, %d)", row, col, row + row_distance, col + col_distance)
    //          << std::endl;
    //   }

    assert(0 <= (nRow + nRowDistance) && (nRow + nRowDistance) <= this->getNumOfRows());
    assert(0 <= (nCol + nColDistance) && (nCol + nColDistance) <= this->getNumOfCols());

    TlMatrix answer(nRowDistance, nColDistance);
#pragma omp parallel for
    for (int dr = 0; dr < nRowDistance; ++dr) {
        const int r = nRow + dr;
        for (int dc = 0; dc < nColDistance; ++dc) {
            const int c = nCol + dc;

            answer.set(dr, dc, this->get(r, c));
        }
    }

    return answer;
}


void TlMatrix::setBlockMatrix(const index_type row, const index_type col, const TlMatrix& matrix)
{
    TlLogging& log = TlLogging::getInstance();
    const int row_distance = matrix.getNumOfRows();
    const int col_distance = matrix.getNumOfCols();

    if (!((0 <= row && row < this->getNumOfRows()) &&
          (0 <= col && col < this->getNumOfCols()) &&
          (0 < (row + row_distance) && (row + row_distance) <= this->getNumOfRows()) &&
          (0 < (col + col_distance) && (col + col_distance) <= this->getNumOfCols()))) {
        log.critical(TlUtils::format("setBlockMatrix() start(%d, %d) mat(%d, %d) -> (%d, %d)",
                                     row, col, matrix.getNumOfRows(), matrix.getNumOfCols(),
                                     this->getNumOfRows(), this->getNumOfCols()));
    }
    assert(0 <= row && row < this->getNumOfRows());
    assert(0 <= col && col < this->getNumOfCols());
    assert(0 < (row + row_distance) && (row + row_distance) <= this->getNumOfRows());
    assert(0 < (col + col_distance) && (col + col_distance) <= this->getNumOfCols());

#pragma omp parallel for
    for (int dr = 0; dr < row_distance; ++dr) {
        const int r = row + dr;
        for (int dc = 0; dc < col_distance; ++dc) {
            const int c = col + dc;

            this->set(r, c,  matrix.get(dr, dc));
        }
    }
}


void TlMatrix::addBlockMatrix(const int row, const int col, const TlMatrix& matrix)
{
    const int row_distance = matrix.getNumOfRows();
    const int col_distance = matrix.getNumOfCols();

    assert(0 <= row && row < this->getNumOfRows());
    assert(0 <= col && col < this->getNumOfCols());
    assert(0 <= (row + row_distance) && (row + row_distance) < this->getNumOfRows());
    assert(0 <= (col + col_distance) && (col + col_distance) < this->getNumOfCols());

#pragma omp parallel for
    for (int dr = 0; dr < row_distance; ++dr) {
        const int r = row + dr;
        for (int dc = 0; dc < col_distance; ++dc) {
            const int c = col + dc;

            this->add(r, c, matrix.get(dr, dc));
        }
    }
}


double TlMatrix::trace() const
{
    assert(this->m_nRows == this->m_nCols);

    const int dim = this->m_nRows;
    double answer = 0.0;
#pragma omp parallel for reduction(+:answer)
    for (int i = 0; i < dim; ++i) {
        answer += this->get(i, i);
    }

    return answer;
}


double TlMatrix::sum() const
{
    //double answer = this->m_aMatrix.sum();
    return std::accumulate(this->data_, this->data_ + this->getNumOfElements(), 0.0);
}


// operator
TlMatrix& TlMatrix::operator=(const TlMatrix& rhs)
{
    if (this != &rhs) {
        this->clear();
        this->m_nRows = rhs.m_nRows;
        this->m_nCols = rhs.m_nCols;
        this->initialize(false);
        const std::size_t size = this->getNumOfElements();
        std::copy(rhs.data_, rhs.data_ + size, this->data_);
    }

    return (*this);
}


TlMatrix& TlMatrix::operator=(const TlSymmetricMatrix& rhs)
{
    this->clear();
    this->m_nRows = rhs.getNumOfRows();
    this->m_nCols = rhs.getNumOfCols();
    this->initialize(false);

    const int numOfRows = this->m_nRows;
    for (int row = 0; row < numOfRows; ++row) {
        // case: row != col
        for (int col = 0; col < row ; ++col) {
            const double tmp = rhs(row, col);
            this->set(row, col, tmp);
            this->set(col, row, tmp);
        }

        // case : row == col
        this->set(row, row, rhs(row, row));
    }

    return (*this);
}


TlMatrix& TlMatrix::operator+=(const TlMatrix& rhs)
{
    assert(this->getNumOfRows() == rhs.getNumOfRows());
    assert(this->getNumOfCols() == rhs.getNumOfCols());

    const size_type size = this->getNumOfElements();

#ifdef _OPENMP
#pragma omp critical(TlMatrix__operator_addplus)
    {
        // use OpenMP
        const size_type quot = size / MAX_LOOP;
        const int rem = size - quot * MAX_LOOP;
#pragma omp parallel
        {
            for (size_type block = 0; block < quot; ++block) {
                const size_type index_base = block * MAX_LOOP;
#pragma omp for
                for (int i = 0; i < MAX_LOOP; ++i) {
                    const size_type index = index_base + i;
                    this->data_[index] += rhs.data_[index];
                }
            }
            
            const size_type index_base = quot * MAX_LOOP;
#pragma omp for
            for (int i = 0; i < rem; ++i) {
                const size_type index = index_base + i;
                this->data_[index] += rhs.data_[index];
            }
        }
    }

#else
    // not use OpenMP
    for (size_type index = 0; index < size; ++index) {
        this->data_[index] += rhs.data_[index];
    }
#endif // _OPENMP

    return (*this);
}


TlMatrix& TlMatrix::operator+=(const TlSymmetricMatrix& rhs)
{
    assert(this->getNumOfRows() == rhs.getNumOfRows());
    assert(this->getNumOfCols() == rhs.getNumOfCols());

    const int numOfRows = rhs.getNumOfRows();

#pragma omp parallel for
    for (int row = 0; row < numOfRows; ++row) {
        // row != col
        for (int col = 0; col < row; ++col) {
            const double tmp = rhs(row, col);
#pragma omp critical(TlMatrix_operator_plus_equal)
            {
                this->add(row, col, tmp);
                this->add(col, row, tmp);
            }
        }

        // row == col
        this->add(row, row, rhs(row, row));
    }

    return (*this);
}


TlMatrix& TlMatrix::operator-=(const TlMatrix& rhs)
{
    assert(this->m_nRows == rhs.m_nRows);
    assert(this->m_nCols == rhs.m_nCols);
    
    const size_type size = this->getNumOfElements();

#ifdef _OPENMP
    // use OpenMP
    const size_type quot = size / MAX_LOOP;
    const int rem = size - quot * MAX_LOOP;
#pragma omp parallel
    {
        for (size_type block = 0; block < quot; ++block) {
            const size_type index_base = block * MAX_LOOP;
#pragma omp for
            for (int i = 0; i < MAX_LOOP; ++i) {
                const size_type index = index_base + i;
                this->data_[index] -= rhs.data_[index];
            }
        }

        const size_type index_base = quot * MAX_LOOP;
#pragma omp for
        for (int i = 0; i < rem; ++i) {
            const size_type index = index_base + i;
            this->data_[index] -= rhs.data_[index];
        }
    }
#else
    for (size_type index = 0; index < size; ++index) {
        this->data_[index] -= rhs.data_[index];
    }
#endif // _OPENMP

    return (*this);
}


TlMatrix& TlMatrix::operator-=(const TlSymmetricMatrix& rhs)
{
    assert(this->getNumOfRows() == rhs.getNumOfRows());
    assert(this->getNumOfCols() == rhs.getNumOfCols());

    const int numOfRows = rhs.getNumOfRows();
#pragma omp parallel for
    for (int row = 0; row < numOfRows; ++row) {
        // row != col
        for (int col = 0; col < row; ++col) {
            const double tmp = - rhs(row, col);
#pragma omp critical(TlMatrix_operator_minus_equal)
            {
                this->add(row, col, tmp);
                this->add(col, row, tmp);
            }
        }

        // row == col
        this->add(row, row, rhs(row, row));
    }

    return (*this);
}


TlMatrix& TlMatrix::operator*=(const TlMatrix& rhs)
{
    TlMatrix tmp = (*this) * rhs;
    (*this) = tmp;

    return (*this);
}


TlMatrix& TlMatrix::operator*=(const double& rhs)
{
    const size_type size = this->getNumOfElements();

#ifdef _OPENMP
    // use OpenMP
    const size_type quot = size / MAX_LOOP;
    const int rem = size - quot * MAX_LOOP;
#pragma omp parallel
    {
        for (size_type block = 0; block < quot; ++block) {
            const size_type index_base = block * MAX_LOOP;
#pragma omp for
            for (int i = 0; i < MAX_LOOP; ++i) {
                const size_type index = index_base + i;
                this->data_[index] *= rhs;
            }
        }

        const size_type index_base = quot * MAX_LOOP;
#pragma omp for
        for (int i = 0; i < rem; ++i) {
            const size_type index = index_base + i;
            this->data_[index] *= rhs;
        }
    }
#else
    // not use OpenMP
    for (size_type index = 0; index < size; ++index) {
        this->data_[index] *= rhs;
    }
#endif // _OPENMP

    return (*this);
}


TlMatrix& TlMatrix::operator/=(const double& rhs)
{
    return (this->operator*=(1.0 / rhs));
}


TlMatrix operator+(const TlMatrix& X, const TlMatrix& Y)
{
    assert(X.m_nRows == Y.m_nRows);
    assert(X.m_nCols == Y.m_nCols);

    TlMatrix answer(X);
    answer += Y;

    return answer;
}


TlMatrix operator-(const TlMatrix& X, const TlMatrix& Y)
{
    assert(X.m_nRows == Y.m_nRows);
    assert(X.m_nCols == Y.m_nCols);

    TlMatrix answer(X);
    answer -= Y;

    return answer;
}


// X x Y
TlMatrix operator*(const TlMatrix& X, const TlMatrix& Y)
{
#ifdef HAVE_LAPACK
    return multiplicationByLapack(X, Y);
#else
    assert(X.m_nCols == Y.m_nRows);
    TlMatrix answer(X.m_nRows, Y.m_nCols);
    for (int row = 0; row < X.m_nRows; ++row) {
        for (int col = 0; col < Y.m_nCols; ++col) {
            for (int t = 0; t < X.m_nCols; ++t) {
                answer(row, col) += X(row, t) * Y(t, col);
            }
        }
    }

    return answer;
#endif // HAVE_LAPACK
}


TlMatrix operator*(const TlMatrix& X, double Y)
{
    TlMatrix answer(X);
    answer *= Y;

    return answer;
}


// 行列x縦ベクトル
TlVector operator*(const TlMatrix& A, const TlVector& X)
{
#ifdef HAVE_LAPACK
    return multiplicationByLapack(A, X);
#else
    {
        assert(A.getNumOfCols() == X.getSize());
        TlVector answer(A.getNumOfRows());
        for (int row = 0; row < X.m_nRows; ++row) {
            double tmp = 0.0;
            for (int col = 0; col < X.m_nCols; ++col) {
                tmp += A.get(row, col) * X.get(col);
            }
            answer.set(row, tmp);
        }
        return answer;
    }
#endif // HAVE_LAPACK
}


// 横ベクトルx行列
TlVector operator*(const TlVector& X, const TlMatrix& A)
{
#ifdef HAVE_LAPACK
    return multiplicationByLapack(X, A);
#else
    {
        assert(A.getNumOfRows() == X.getSize());
        TlVector answer(A.getNumOfCol());
        for (int col = 0; col < X.m_nCols; ++col) {
            double tmp = 0.0;
            for (int row = 0; row < X.m_nRows; ++row) {
                tmp += A.get(row, col) * X.get(row);
            }
            answer.set(col, tmp);
        }
        return answer;
    }
#endif // HAVE_LAPACK
}

// ddot?
const TlMatrix& TlMatrix::dot(const TlMatrix& X)
{
    assert(this->getNumOfRows() == X.getNumOfRows());
    assert(this->getNumOfCols() == X.getNumOfCols());

    const size_type size = this->getNumOfElements();

#ifdef _OPENMP
    // use OpenMP
    const size_type quot = size / MAX_LOOP;
    const int rem = size - quot * MAX_LOOP;
#pragma omp parallel
    {
        for (size_type block = 0; block < quot; ++block) {
            const size_type index_base = block * MAX_LOOP;
#pragma omp for
            for (int i = 0; i < MAX_LOOP; ++i) {
                const size_type index = index_base + i;
                this->data_[index] *= X.data_[index];
            }
        }

        const size_type index_base = quot * MAX_LOOP;
#pragma omp for
        for (int i = 0; i < rem; ++i) {
            const size_type index = index_base + i;
            this->data_[index] *= X.data_[index];
        }
    }
#else
    // not use OpenMP
    for (size_type index = 0; index < size; ++index) {
        this->data_[index] *= X.data_[index];
    }
#endif // _OPENMP
    
    return (*this);
}


const TlMatrix& TlMatrix::dot(const TlSymmetricMatrix& X)
{
    this->dot(TlMatrix(X));

    return *this;
}


TlMatrix dot(const TlMatrix& X, const TlMatrix& Y)
{
    assert(X.getNumOfRows() == Y.getNumOfRows());
    assert(X.getNumOfCols() == Y.getNumOfCols());

    TlMatrix Z = X;
    Z.dot(Y);

    return Z;
}


std::size_t TlMatrix::getMemSize() const
{
    std::size_t answer = sizeof(TlMatrix);
    answer += sizeof(double) * this->getNumOfElements();

    return answer;
}   


/**
   入力されたstd::ifstream が読み込み可能なフォーマットかチェックする。
   チェック後、渡されたifstream の読み込み位置は先頭に戻される。
 */
bool TlMatrix::isLoadable(std::ifstream& ifs)
{
    if (ifs.is_open() != true) {
        return false;
    }

    enum {RSFD, CSFD, RLHD, CLHD, RUHD, CUHD, RSFS, CSFS, RLHS, CLHS, RUHS, CUHS, LHSC};
    int nType = 0;

    ifs.seekg(0, std::ios::beg);
    const bool bHeader = TlMatrix::getHeaderInfo(ifs, &nType);
    ifs.seekg(0, std::ios::beg);

    bool bJudge = false;
    if (bHeader == true) {
        if (nType == RSFD) {
            bJudge =  true;
        }
    }

    return bJudge;
}


bool TlMatrix::isLoadable(const std::string& rFilePath)
{
    std::ifstream ifs;
    ifs.open(rFilePath.c_str());
    if (ifs.fail()) {
#ifdef DEBUG
        std::cerr << "[error] TlMatrix::load(): could not open file. " << rFilePath << std::endl;
#endif // DEBUG
        return false;
    }

    const bool bAnswer = TlMatrix::isLoadable(ifs);
    ifs.close();
    return bAnswer;
}


bool TlMatrix::getHeaderInfo(std::ifstream& ifs, int* pType,
                             int* pNumOfRows, int* pNumOfCols)
{
    bool bAnswer = true;

    // get file size
    std::ifstream::pos_type nFileSize = 0;
    {
        ifs.seekg(0, std::ios_base::beg);
        std::ifstream::pos_type begin = ifs.tellg();
        ifs.seekg(0, std::ios_base::end);
        std::ifstream::pos_type end = ifs.tellg();

        nFileSize = end - begin;
    }

    int  nType = 0;
    bool bCheckVariableType = false;
    std::ifstream::pos_type nStartContentPos = 0;

    // int case:
    {
        int nRows = 0;
        int nCols = 0;
        ifs.seekg(0, std::ios_base::beg);
        ifs.read((char*)&(nType), sizeof(int));
        ifs.read((char*)&(nRows), sizeof(int));
        ifs.read((char*)&(nCols), sizeof(int));

        const std::ifstream::pos_type nEstimatedFileSize =
            std::ifstream::pos_type(sizeof(int) *3)
            + (std::ifstream::pos_type(nRows) * std::ifstream::pos_type(nCols) * std::ifstream::pos_type(sizeof(double)));
        if (nEstimatedFileSize == nFileSize) {
            if (pType != NULL) {
                *pType = nType;
            }
            if (pNumOfRows != NULL) {
                *pNumOfRows = nRows;
            }
            if (pNumOfCols != NULL) {
                *pNumOfCols = nCols;
            }
            bCheckVariableType = true;
            nStartContentPos = sizeof(int) * 3;
        }
    }

    // long case:
    {
        long nRows = 0;
        long nCols = 0;
        ifs.seekg(0, std::ios_base::beg);
        ifs.read((char*)&(nType), sizeof(int));
        ifs.read((char*)&(nRows), sizeof(long));
        ifs.read((char*)&(nCols), sizeof(long));

        const std::ifstream::pos_type nEstimatedFileSize =
            std::ifstream::pos_type(sizeof(int) + sizeof(long) * 2)
            + (std::ifstream::pos_type(nRows) * std::ifstream::pos_type(nCols) * std::ifstream::pos_type(sizeof(double)));
        if (nEstimatedFileSize == nFileSize) {
            if (pType != NULL) {
                *pType = nType;
            }
            if (pNumOfRows != NULL) {
                *pNumOfRows = static_cast<int>(nRows);
            }
            if (pNumOfCols != NULL) {
                *pNumOfCols = static_cast<int>(nCols);
            }
            bCheckVariableType = true;
            nStartContentPos = sizeof(int) + sizeof(long) * 2;
        }
    }

    if (bCheckVariableType == true) {
        ifs.seekg(nStartContentPos, std::ios_base::beg);
    } else {
        //std::cerr << "file size mismatch." << std::endl;
        bAnswer = false;
    }

    return bAnswer;
}


bool TlMatrix::load(const std::string& sFilePath)
{
    std::ifstream ifs;

    ifs.open(sFilePath.c_str());
    if (ifs.fail()) {
#ifdef DEBUG
        std::cerr << "[error] TlMatrix::load(): could not open file. " << sFilePath << std::endl;
#endif // DEBUG
        return false;
    }

    bool bAnswer = this->load(ifs);
    ifs.close();

    if (bAnswer != true) {
        std::cerr << "TlMatrix::load() is not supported: " << sFilePath << std::endl;
        return false;
    }

    return true;
}


bool TlMatrix::load(std::ifstream& ifs)
{
    bool bAnswer = true;

    // binary mode
    {
        enum {RSFD, CSFD, RLHD, CLHD, RUHD, CUHD, RSFS, CSFS, RLHS, CLHS, RUHS, CUHS, LHSC};

        // read header
        int nType = 0;
        int nRows = 0;
        int nCols = 0;

        bAnswer = TlMatrix::getHeaderInfo(ifs, &nType, &nRows, &nCols);
        switch (nType) {
        case RSFD:
            //     std::cerr << "load RSFD" << std::endl;
            break;

        default:
            std::cerr << "this matrix type is not supported." << std::endl;
            bAnswer = false;
            break;
        }

        if (bAnswer == true) {
            this->clear();
            this->m_nRows = nRows;
            this->m_nCols = nCols;
            this->initialize();

            // 本当は大きめのバッファで読み取りたいが、row oriented だから変換しないといけない
            //const int nMatrixSize = this->m_nRows * this->m_nCols;
            //ifs.read(reinterpret_cast<char*>(&(this->m_aMatrix[0])), sizeof(double) * nMatrixSize);
            for (int row = 0; row < nRows; ++row) {
                for (int col = 0; col < nCols; ++col) {
                    double tmp;
                    ifs.read((char*)&tmp, sizeof(double));
                    (*this)(row, col) = tmp;
                }
            }
        }
    }

    return bAnswer;
}


bool TlMatrix::save(const std::string& sFilePath) const
{
    bool bAnswer = true;

    std::ofstream ofs;
    ofs.open(sFilePath.c_str(), std::ofstream::out | std::ofstream::binary);

    bAnswer = this->save(ofs);

    ofs.close();

    return bAnswer;
}


bool TlMatrix::save(std::ofstream& ofs) const
{
    bool bAnswer = true;

    const int nType = 0; // means RSFD
    ofs.write(reinterpret_cast<const char*>(&nType), sizeof(int));
    ofs.write(reinterpret_cast<const char*>(&this->m_nRows), sizeof(int));
    ofs.write(reinterpret_cast<const char*>(&this->m_nCols), sizeof(int));

    // 本当は大きめのバッファで書き込みたいが、RSFD(row oriented)だから変換しないといけない
    const index_type maxRow = this->m_nRows;
    const index_type maxCol = this->m_nCols;
    for (index_type row = 0; row < maxRow; ++row) {
        for (index_type col = 0; col < maxCol; ++col) {
            const double tmp = this->get(row, col);
            ofs.write(reinterpret_cast<const char*>(&tmp), sizeof(double));
        }
    }

    return bAnswer;
}


bool TlMatrix::saveText(const std::string& sFilePath) const
{
    bool bAnswer = true;

    std::ofstream ofs;
    ofs.open(sFilePath.c_str(), std::ofstream::out);

    if (ofs.good()) {
        bAnswer = this->saveText(ofs);
    } else {
        bAnswer =false;
    }

    ofs.close();

    return bAnswer;
}


bool TlMatrix::saveText(std::ostream& os) const
{
    bool bAnswer = true;
    const int nRows = this->getNumOfRows();
    const int nCols = this->getNumOfCols();

    os << "TEXT\n";
    os << nRows << "\n";
    os << nCols << "\n";

    // print out LCAO coefficent
    for (int i=0; i < nRows; ++i) {
        for (int j=0; j < nCols; ++j) {
            os << TlUtils::format(" %10.6lf", (*this)(i, j));
        }
        os << "\n";
    }
    os << "\n";

    return bAnswer;
}


TlSerializeData TlMatrix::getSerialize() const
{
    TlSerializeData data;
    data["row"] = this->getNumOfRows();
    data["col"] = this->getNumOfCols();
    data["type"] = "RSFD";

    TlSerializeData tmp;
    const size_type size = this->getNumOfElements();
    for (size_type index = 0; index < size; ++index) {
        tmp.pushBack(this->data_[index]);
    }
    data["data"] = tmp;

    return data;
}


/**
 *  出力ストリーム
 */
std::ostream& operator <<(std::ostream& out, const TlMatrix& rhs)
{
    rhs.print(out);
    return out;
}


/**
 *  CSV出力
 */
std::string TlMatrix::getCsv() const
{
    std::ostringstream out;

    for (int row = 0; row < this->m_nRows; ++row) {
        for (int col = 0; col < this->m_nCols; ++col) {
            out << (*this)(row, col) << ", ";
        }
        out << std::endl;
    }

    return out.str();
}

// 転置
// TODO: 高速化できる？
const TlMatrix& TlMatrix::transpose()
{
    const int numOfRows = this->getNumOfRows();
    const int numOfCols = this->getNumOfCols();

    if (numOfRows == numOfCols) {
        // 正方行列の場合
        for (int row = 0; row < numOfRows; ++row) {
            for (int col = 0; col < row; ++col) {
                std::swap<double>((*this)(row, col), (*this)(col, row));
            }
        }
    } else {
        // 正方行列でない場合
        TlMatrix tmp(*this);
        this->resize(numOfCols, numOfRows);
        for (int row = 0; row < numOfRows; ++row) {
            for (int col = 0; col < numOfCols; ++col) {
                (*this)(col, row) = tmp(row, col);
            }
        }
    }

    return (*this);
}


bool TlMatrix::inverse()
{
#ifdef HAVE_LAPACK
    // using LAPACK
    return inverseByLapack(*this);
#else
    // without LAPACK
    std::cerr << "sorry. this code is not implemented." << std::endl;
    abort();
    return false;
#endif // HAVE_LAPACK  
}


TlVector TlMatrix::getDiagonalElements() const
{
    const index_type dim = std::min(this->getNumOfRows(), this->getNumOfCols());
    TlVector answer(dim);
    for (index_type i = 0; i < dim; ++i) {
        const double value = this->get(i, i);
        answer.set(i, value);
    }

    return answer;
}



// 要素の絶対値の最大値を返す
// int* outRow, outCol にはその要素がある行と列を返す
double TlMatrix::getMaxAbsoluteElement(int* outRow, int* outCol) const
{
    double dAnswer = -1.0;
    int nMaxRow = 0;
    int nMaxCol = 0;

    const int numOfRows = this->getNumOfRows();
    const int numOfCols = this->getNumOfCols();
    for (int row = 0; row < numOfRows; ++row) {
        for (int col = 0; col < numOfCols; ++col) {
            double val = fabs((*this)(row, col));
            if (dAnswer < val) {
                dAnswer = val;
                nMaxRow = row;
                nMaxCol = col;
            }
        }
    }

    if (outRow != NULL) {
        *outRow = nMaxRow;
    }
    if (outCol != NULL) {
        *outCol = nMaxCol;
    }

    return dAnswer;
}


#ifdef HAVE_LAPACK
extern "C" {
    void dgemv_(const char* TRANS, const int* M, const int* N,
                const double* ALPHA, const double* A, const int* LDA,
                const double* X, const int* INCX, const double* BETA,
                double* Y, const int* INCY);
    
    void dgemm_(const char*, const char*, const int*, const int*, const int*,
    const double*, const double*, int*, const double*, const int* ,const double*, double*, int*);

    void dgetrf_(const int* M, const int* N, double* A, const int* LDA, int* IPIV, int* INFO);
    void dgetri_(const int* N, double* A, const int* LDA, int* IPIV, double* WORK, int* LWORK, int* INFO);
}


TlVector multiplicationByLapack(const TlMatrix& A, const TlVector& X)
{
    if (A.getNumOfCols() != X.getSize()) {
        std::cerr << TlUtils::format("multiplicationByLapack(): parameter error at file:%s, line:%d",
                                     __FILE__, __LINE__)
                  << std::endl;
        std::abort();
    }

    TlVector Y(A.getNumOfRows());
    const char TRANS = 'N';
    const int M = A.getNumOfRows();
    const int N = A.getNumOfCols();
    const double alpha = 1.0;
    const int LDA = std::max(1, M);
    const int INCX = 1;
    const double beta = 1.0;
    const int INCY = 1;
    
    // *  DGEMV  performs one of the matrix-vector operations
    // *
    // *     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
    // *
    // *  where alpha and beta are scalars, x and y are vectors and A is an
    // *  m by n matrix.
    dgemv_(&TRANS, &M, &N,
           &alpha, const_cast<TlMatrix&>(A).data_, &LDA,
           const_cast<TlVector&>(X).data_, &INCX,
           &beta, Y.data_, &INCY);

    return Y;
}


TlVector multiplicationByLapack(const TlVector& X, const TlMatrix& A)
{
    if (X.getSize() != A.getNumOfRows()) {
        std::cerr << TlUtils::format("multiplicationByLapack(): parameter error at file:%s, line:%d",
                                     __FILE__, __LINE__)
                  << std::endl;
        std::abort();
    }

    TlVector Y(A.getNumOfCols());
    const char TRANS = 'T';
    const int M = A.getNumOfRows();
    const int N = A.getNumOfCols();
    const double alpha = 1.0;
    const int LDA = std::max(1, M);
    const int INCX = 1;
    const double beta = 1.0;
    const int INCY = 1;
    
    // *  DGEMV  performs one of the matrix-vector operations
    // *
    // *     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
    // *
    // *  where alpha and beta are scalars, x and y are vectors and A is an
    // *  m by n matrix.
    dgemv_(&TRANS, &M, &N,
           &alpha, const_cast<TlMatrix&>(A).data_, &LDA,
           const_cast<TlVector&>(X).data_, &INCX,
           &beta, Y.data_, &INCY);

    return Y;
}


TlMatrix multiplicationByLapack(const TlMatrix& X, const TlMatrix& Y)
{
    //std::cerr << "use dgemm() of LAPACK" << std::endl;
    if (X.getNumOfCols() != Y.getNumOfRows()) {
        std::cerr << TlUtils::format("multiplicationByLapack(): parameter error at file:%s, line:%d",
                                     __FILE__, __LINE__)
                  << std::endl;
        std::abort();
    }

    TlMatrix answer(X.m_nRows, Y.m_nCols);

    char TRANSA = 'N';
    char TRANSB = 'N';
    double alpha = 1.0;
    double beta  = 0.0;
    int M = X.m_nRows;
    int N = Y.m_nCols;
    int K = X.m_nCols;
    int LDA = M;
    int LDB = K;
    int LDC = M;

    // DGEMM  performs one of the matrix-matrix operations
    // C := alpha*op( A )*op( B ) + beta*C,
    // where  op( X ) is one of op( X ) = X   or   op( X ) = X',
    // alpha and beta are scalars, and A, B and C are matrices, with op( A )
    // an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix
    dgemm_(&TRANSA, &TRANSB, &M, &N, &K, &alpha,
           const_cast<TlMatrix&>(X).data_, &LDA,
           const_cast<TlMatrix&>(Y).data_, &LDB, &beta,
           answer.data_, &LDC);

    return answer;
}

bool inverseByLapack(TlMatrix& X)
{
    bool bAnswer = false;
    
    const int M = X.getNumOfRows();
    const int N = X.getNumOfCols();

    const int LDA = std::max(1, M);
    int* IPIV = new int[std::min(M, N)];
    int INFO = 0;

    double* A = X.data_;
    int LWORK = std::max(1, N);
    double* WORK = new double[LWORK];
    
    dgetrf_(&M, &N, A, &LDA, IPIV, &INFO);
    if (INFO == 0) {
        dgetri_(&N, A, &LDA, IPIV, WORK, &LWORK, &INFO);
        if (INFO == 0) {
            bAnswer = true;
        } else {
            std::cerr << "inverseByLapack() failed.: dgetri() return code = " << INFO << std::endl;
        }
    } else {
        std::cerr << "inverseByLapack() failed.: dgetrf() return code = " << INFO << std::endl;
    }

    std::swap(X.m_nRows, X.m_nCols);
    
    delete[] WORK;
    WORK = NULL;
    
    delete[] IPIV;
    IPIV = NULL;

    return bAnswer;
}

#endif // HAVE_LAPACK

