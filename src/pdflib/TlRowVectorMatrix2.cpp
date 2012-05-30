#include <iostream>
#include "TlRowVectorMatrix2.h"
#include "TlCommunicate.h"

TlRowVectorMatrix2::TlRowVectorMatrix2(const index_type row,
                                       const index_type col,
                                       int allProcs, int rank)
    : numOfRows_(0), numOfCols_(0), reserveCols_(0),
      allProcs_(allProcs), rank_(rank),
      numOfLocalRows_(0) {
    this->resize(row, col);
}


TlRowVectorMatrix2::~TlRowVectorMatrix2()
{
}


void TlRowVectorMatrix2::resize(const index_type newRows,
                                const index_type newCols)
{
    this->numOfRows_ = newRows;
    this->numOfCols_ = newCols;

    const div_t turns = std::div(newRows, this->allProcs_);
    index_type numOfLocalRows = turns.quot;
    if (this->rank_ < turns.rem) {
        numOfLocalRows += 1;
    }
    this->numOfLocalRows_ = numOfLocalRows;
    this->data_.resize(numOfLocalRows);

    if (this->reserveCols_ < newCols) {
        this->reserveCols_ = newCols;
    }
    for (index_type i = 0; i < numOfLocalRows; ++i) {
        this->data_[i].reserve(this->reserveCols_);
        this->data_[i].resize(newCols);
    }
}


void TlRowVectorMatrix2::reserve_cols(const index_type newReserves) {
    this->reserveCols_ = std::max(this->getNumOfCols(), newReserves);
}


void TlRowVectorMatrix2::set(index_type row, index_type col, double value)
{
    const div_t turns = std::div(row, this->allProcs_);
    if (turns.rem == this->rank_) {
        const index_type row = turns.quot;
        this->data_[row][col] = value;
    }
}


TlVector TlRowVectorMatrix2::getRowVector(index_type globalRow) const
{
    TlVector answer;
    const div_t turns = std::div(globalRow, this->allProcs_);
    if (turns.rem == this->rank_) {
        const index_type row = turns.quot;
        // std::cerr << globalRow << ", " 
        //           << row << ", " << this->numOfLocalRows_ << ", " 
        //           << this->data_.size() << std::endl;
        assert(row < this->numOfLocalRows_);
        answer = TlVector(this->data_[row]);
    }

    return answer;
}


TlMatrixObject::index_type 
TlRowVectorMatrix2::getRowVector(index_type row,
                                double *pBuf,
                                index_type maxColSize) const
{
    index_type copySize = 0;
    const div_t turns = std::div(row, this->allProcs_);
    if (turns.rem == this->rank_) {
        const index_type row = turns.quot;
        assert(row < this->numOfLocalRows_);

        copySize = std::min(this->getNumOfCols(), maxColSize);
        std::copy(this->data_[row].begin(),
                  this->data_[row].begin() + copySize,
                  pBuf);
    }

    return copySize;
}


int TlRowVectorMatrix2::getPEinChargeByRow(const index_type row) const
{
    assert((0 <= row) && (row < this->getNumOfRows()));
    const div_t turns = std::div(row, this->allProcs_);
    return turns.rem;
}


TlMatrix TlRowVectorMatrix2::getTlMatrix() const 
{
    const index_type numOfRows = this->getNumOfRows();
    const index_type numOfCols = this->getNumOfCols();
    TlMatrix answer(numOfRows, numOfCols);

    // const div_t turns = std::div(numOfRows, this->allProcs_);
    // const index_type localRows = turns.quot + ((this->rank_ < turns.rem) ? 1 : 0);
    const index_type numOfLocalRows = this->numOfLocalRows_;
    for (index_type r = 0; r < numOfLocalRows; ++r) {
        const index_type row = r * this->allProcs_ + this->rank_;
        for (index_type col = 0; col < numOfCols; ++col) {
            answer.set(row, col, this->data_[r][col]);
        }
    }
    
    return answer;
}




