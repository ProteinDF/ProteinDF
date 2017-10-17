// Copyright (C) 2002-2014 The ProteinDF project
// see also AUTHORS and README.
// 
// This file is part of ProteinDF.
// 
// ProteinDF is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// ProteinDF is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

#include "TlSparseVectorMatrix.h"

#define MINIMUM_CACHING_INDECES (2) // 最低保持するキャッシュ数

TlSparseVectorMatrix::TlSparseVectorMatrix(index_type row, index_type col,
                                           std::size_t cacheMemSize)
    : TlMatrixObject(RSFD), numOfRows_(row), numOfCols_(col), cacheMemSize_(cacheMemSize) {
    this->data_.clear();
}


TlSparseVectorMatrix::TlSparseVectorMatrix(const TlSparseVectorMatrix& rhs)
    : TlMatrixObject(RSFD), numOfRows_(rhs.numOfRows_), numOfCols_(rhs.numOfCols_), data_(rhs.data_) {
}


TlSparseVectorMatrix::~TlSparseVectorMatrix()
{
}


std::size_t TlSparseVectorMatrix::getMemSize() const
{
    std::size_t answer = sizeof(index_type) * 2;
    answer += this->data_.size() * sizeof(double) * this->getNumOfCols();

    return answer;
}


TlMatrixObject::index_type TlSparseVectorMatrix::getNumOfRows() const
{
    return this->numOfRows_;
}


TlMatrixObject::index_type TlSparseVectorMatrix::getNumOfCols() const
{
    return this->numOfCols_;
}


double TlSparseVectorMatrix::get(index_type row, index_type col) const
{
    double answer = 0.0;
    DataType::const_iterator it = this->data_.find(row);
    if (it != this->data_.end()) {
        answer = it->second.get(col);
    }

    return answer;
}


void TlSparseVectorMatrix::set(index_type row, index_type col, double value)
{
    DataType::iterator it = this->data_.lower_bound(row);
    if ((it == this->data_.end()) || (it->first != row)) {
        TlVector tmp(this->getNumOfCols());
        this->data_.insert(it, std::pair<index_type, TlVector>(row, tmp));
        it = this->data_.find(row);
    }

    it->second.set(col, value);
}


void TlSparseVectorMatrix::add(index_type row, index_type col, double value)
{
    DataType::iterator it = this->data_.lower_bound(row);
    if ((it == this->data_.end()) || (it->first != row)) {
        TlVector tmp(this->getNumOfCols());
        this->data_.insert(it, std::pair<index_type, TlVector>(row, tmp));
        it = this->data_.find(row);
    }

    it->second.add(col, value);
}


TlVector TlSparseVectorMatrix::getRowVector(index_type row) const
{
    TlVector answer(this->getNumOfCols());
    DataType::const_iterator it = this->data_.find(row);
    if (it != this->data_.end()) {
        answer = it->second;
    }
    
    return answer;
}


TlVector TlSparseVectorMatrix::getColVector(index_type col) const
{
    const index_type numOfRows = this->getNumOfRows();
    TlVector answer(numOfRows);
    
    for (index_type r = 0; r < numOfRows; ++r) {
        answer[r] = this->get(r, col);
    }

    return answer;
}


void TlSparseVectorMatrix::setRowVector(const index_type row, const TlVector& cols,
                                        bool isKeeped)
{
    this->data_[row] = cols;

    if (isKeeped == false) {
        this->updateCache(row);
    }
}


bool TlSparseVectorMatrix::findRow(index_type row) const
{
    DataType::const_iterator it = this->data_.find(row);
    return (it != this->data_.end());
}


void TlSparseVectorMatrix::updateCache(const std::size_t row)
{
    // update cacheIndexList_
    {
        std::list<int>::iterator it = std::find(this->cacheIndexList_.begin(),
                                                this->cacheIndexList_.end(),
                                                row);
        if (it != this->cacheIndexList_.end()) {
            this->cacheIndexList_.erase(it);
        }
        this->cacheIndexList_.push_back(row);
    }

    // remove top index of cacheIndexList_
    const int numOfCachingIndeces = this->cacheIndexList_.size();
    if (numOfCachingIndeces > MINIMUM_CACHING_INDECES) {
        const std::size_t currentMemSize = this->data_.size() * this->getNumOfCols() * sizeof(double);
        if (currentMemSize > this->cacheMemSize_) {
            const int removeRow = this->cacheIndexList_.front();
            
            this->cacheIndexList_.erase(this->cacheIndexList_.begin());
            DataType::iterator it = this->data_.find(removeRow);
            assert(it != this->data_.end());
            this->data_.erase(it);
        }
    }
}
