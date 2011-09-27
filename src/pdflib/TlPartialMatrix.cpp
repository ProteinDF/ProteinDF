#include <iostream>
#include <cassert>
#include "TlPartialMatrix.h"

TlPartialMatrix::TlPartialMatrix(int globalNumOfRows, int globalNumOfCols,
                                 int startRow, int startCol,
                                 int rowRange, int colRange)
    : numOfRows_(globalNumOfRows), numOfCols_(globalNumOfCols),
      startRow_(startRow), startCol_(startCol), pData_(NULL)
{
    this->localNumOfRows_ = std::min(rowRange, this->numOfRows_ - startRow);
    this->localNumOfCols_ = std::min(colRange, this->numOfCols_ - startCol);

    const std::size_t bufCount = this->getRowRange() * this->getColRange();
    this->pData_ = new double[bufCount];
    for (std::size_t i = 0; i < bufCount; ++i) {
        this->pData_[i] = 0.0;
    }
}


TlPartialMatrix::TlPartialMatrix(const TlPartialMatrix& rhs)
    : numOfRows_(rhs.numOfRows_), numOfCols_(rhs.numOfCols_),
      startRow_(rhs.startRow_), startCol_(rhs.startCol_),
      localNumOfRows_(rhs.localNumOfRows_), localNumOfCols_(rhs.localNumOfCols_), pData_(NULL)
{
    const std::size_t bufCount = this->getRowRange() * this->getColRange();
    this->pData_ = new double[bufCount];
    std::copy(rhs.pData_, rhs.pData_ + bufCount, this->pData_);
}


TlPartialMatrix::~TlPartialMatrix()
{
    if (this->pData_ != NULL) {
        delete[] this->pData_;
        this->pData_ = NULL;
    }
}


TlPartialMatrix& TlPartialMatrix::operator=(const TlPartialMatrix& rhs)
{
    if (&rhs != this) {
        this->numOfRows_ = rhs.numOfRows_;
        this->numOfCols_ = rhs.numOfCols_;
        this->startRow_ = rhs.startRow_;
        this->startCol_ = rhs.startCol_;
        this->localNumOfRows_ = rhs.localNumOfRows_;
        this->localNumOfCols_ = rhs.localNumOfCols_;

        if (this->pData_ != NULL) {
            delete[] this->pData_;
            this->pData_ = NULL;
        }
        const std::size_t bufSize = this->getRowRange() * this->getColRange();
        this->pData_ = new double[bufSize];
        std::copy(rhs.pData_, rhs.pData_ + bufSize, this->pData_);
    }

    return *this;
}


TlMatrixObject::index_type TlPartialMatrix::getNumOfRows() const
{
    return this->numOfRows_;
}


TlMatrixObject::index_type TlPartialMatrix::getNumOfCols() const
{
    return this->numOfCols_;
}


TlMatrixObject::index_type TlPartialMatrix::getStartRow() const
{
    return this->startRow_;
}


TlMatrixObject::index_type TlPartialMatrix::getStartCol() const
{
    return this->startCol_;
}


std::size_t TlPartialMatrix::getMemSize() const
{
    std::size_t answer = sizeof(TlPartialMatrix);

    const std::size_t bufSize = this->getRowRange() * this->getColRange();
    answer += sizeof(double) * bufSize ;

    return answer;
}


void TlPartialMatrix::set(index_type globalRow, index_type globalCol, const double value)
{
    const size_type index = this->index(globalRow, globalCol);
    if (index != -1) {
        this->pData_[index] = value;
    }
}


void TlPartialMatrix::add(index_type globalRow, index_type globalCol, const double value)
{
    const size_type index = this->index(globalRow, globalCol);
    if (index != -1) {
        this->pData_[index] += value;
    }
}


double TlPartialMatrix::get(const index_type globalRow, const index_type globalCol) const
{
    double answer = 0.0;
    const size_type index = this->index(globalRow, globalCol);
    if (index != -1) {
        answer = this->pData_[index];
    }

    return answer;
}


double TlPartialMatrix::getLocal(const index_type localRow, const index_type localCol) const
{
    assert(0 <= localRow);
    assert(localRow < this->localNumOfRows_);
    assert(0 <= localCol);
    assert(localCol < this->localNumOfCols_);

    return this->pData_[localRow * this->localNumOfCols_ + localCol];
}


double TlPartialMatrix::getMaxAbsoluteElement(index_type* outRow, index_type* outCol) const
{
    const int max_r = this->getRowRange();
    const int max_c = this->getColRange();

    int row = 0;
    int col = 0;
    double answer = 0.0;
    for (int r = 0; r < max_r; ++r) {
        for (int c = 0; c < max_c; ++c) {
            double tmp = this->getLocal(r, c);
            if (answer < tmp) {
                answer = tmp;
                row = r;
                col = c;
            }
        }
    }

    if (outRow != NULL) {
        *outRow = row + this->startRow_;
    }
    if (outCol != NULL) {
        *outCol = col + this->startCol_;
    }

    return answer;
}


TlVector TlPartialMatrix::getRowVector(const index_type row) const
{
    assert((0 <= row) && (row < this->getNumOfRows()));

    const index_type numOfCols = this->getNumOfCols();
    TlVector answer(numOfCols);

    const index_type startRow = this->getStartRow();
    const index_type endRow = startRow + this->getRowRange();
    const index_type startCol = this->getStartCol();
    const index_type endCol = startCol + this->getColRange();
    if ((startRow <= row) && (row < endRow)) {
#pragma omp parallel for
        for (index_type col = startCol; col < endCol; ++col) {
            answer[col] = this->getLocal(row, col);
        }
    }
        
    return answer;
}


TlVector TlPartialMatrix::getColVector(const index_type col) const
{
    assert((0 <= col) && (col < this->getNumOfCols()));

    const index_type numOfCols = this->getNumOfRows();
    TlVector answer(numOfCols);

    const index_type startRow = this->getStartRow();
    const index_type endRow = startRow + this->getRowRange();
    const index_type startCol = this->getStartCol();
    const index_type endCol = startCol + this->getColRange();
    if ((startCol <= col) && (col < endCol)) {
#pragma omp parallel for
        for (index_type row = startRow; row < endRow; ++row) {
            answer[row] = this->getLocal(row, col);
        }
    }
        
    return answer;
}
