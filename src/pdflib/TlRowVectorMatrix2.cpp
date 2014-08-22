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

#include <iostream>
#include <new>
#include "TlRowVectorMatrix2.h"
#include "TlMemManager.h"

TlRowVectorMatrix2::TlRowVectorMatrix2(const index_type row,
                                       const index_type col,
                                       int allProcs, int rank,
                                       bool isUsingMemManager)
    : log_(TlLogging::getInstance()),
      numOfRows_(0), numOfCols_(0), reserveCols_(0),
      allProcs_(allProcs), rank_(rank),
      numOfLocalRows_(0), isUsingMemManager_(isUsingMemManager) {
    this->resize(row, col);
}

TlRowVectorMatrix2::TlRowVectorMatrix2(const TlRowVectorMatrix2& rhs)
    : log_(TlLogging::getInstance()),
      numOfRows_(0), numOfCols_(0), reserveCols_(0),
      allProcs_(rhs.allProcs_), rank_(rhs.rank_),
      numOfLocalRows_(0), isUsingMemManager_(rhs.isUsingMemManager_) {

    this->resize(rhs.getNumOfRows(), rhs.getNumOfCols());

    assert(this->data_.size() == rhs.data_.size());
    const index_type numOfCols = this->getNumOfCols();
    const std::size_t size = this->data_.size();
    for (std::size_t i = 0; i < size; ++i) {
        assert(this->data_[i] != NULL);
        std::copy(rhs.data_[i], rhs.data_[i] + numOfCols, this->data_[i]);
    }
}


TlRowVectorMatrix2::~TlRowVectorMatrix2()
{
    const index_type numOfLocalRows = this->numOfLocalRows_;

    if (this->isUsingMemManager_ == true) {
        TlMemManager& rMemManager = TlMemManager::getInstance();
        const index_type reserveCols = this->reserveCols_;
        for (index_type i = 0; i < numOfLocalRows; ++i) {
            rMemManager.deallocate((char*)this->data_[i]);
            this->data_[i] = NULL;
        }
    } else {
        for (index_type i = 0; i < numOfLocalRows; ++i) {
            delete[] this->data_[i];
            this->data_[i] = NULL;
        }
    }
    this->data_.clear();
}


void TlRowVectorMatrix2::resize(const index_type newRows,
                                const index_type newCols)
{
    this->numOfRows_ = newRows;
    index_type prevNumOfLocalRows = this->numOfLocalRows_;
    const div_t turns = std::div(newRows, this->allProcs_);
    index_type newNumOfLocalRows = turns.quot;
    if (this->rank_ < turns.rem) {
        newNumOfLocalRows += 1;
    }
    this->numOfLocalRows_ = newNumOfLocalRows;
    if (newNumOfLocalRows > prevNumOfLocalRows) {
        this->data_.resize(newNumOfLocalRows, NULL);
    } else if (newNumOfLocalRows < prevNumOfLocalRows) {
        if (this->isUsingMemManager_ == true) {
            TlMemManager& rMemManager = TlMemManager::getInstance();
            const index_type reserveCols = this->reserveCols_;
            for (index_type i = newNumOfLocalRows; i < prevNumOfLocalRows; ++i) {
                rMemManager.deallocate((char*)this->data_[i]);
                this->data_[i] = NULL;
            }
        } else {
            for (index_type i = newNumOfLocalRows; i < prevNumOfLocalRows; ++i) {
                delete[] this->data_[i];
                this->data_[i] = NULL;
            }
        }
        this->data_.resize(newNumOfLocalRows);
    }
    
    this->numOfCols_ = newCols;
    this->reserve_cols(newCols);
}


void TlRowVectorMatrix2::reserve_cols(const index_type newReserves) {
    const index_type prevReserveCols = this->reserveCols_;
    const index_type newReserveCols = std::max(this->getNumOfCols(), newReserves);
    // const index_type numOfCols = this->getNumOfCols();

    if (prevReserveCols < newReserveCols) {
        this->reserveCols_ = newReserveCols;
        const index_type numOfLocalRows = this->numOfLocalRows_;
        for (index_type i = 0; i < numOfLocalRows; ++i) {
            double* pNew = NULL;
            try {
                if (this->isUsingMemManager_ == true) {
                    TlMemManager& rMemManager = TlMemManager::getInstance();
                    pNew = (double*)rMemManager.allocate(sizeof(double)*newReserveCols);
                } else {
                    pNew = new double[newReserveCols];
                }
            } catch (std::bad_alloc& ba) {
                this->log_.critical(TlUtils::format("bad_alloc caught: %s", ba.what()));
                throw;
            } catch (...) {
                this->log_.critical("unknown error.");
                throw;
            }
            assert(pNew != NULL);

            for (index_type j = 0; j < newReserveCols; ++j) {
                pNew[j] = 0.0;
            }
            
            if (this->data_[i] != NULL) {
                for (index_type j = 0; j < prevReserveCols; ++j) {
                    pNew[j] = this->data_[i][j];
                }
                if (this->isUsingMemManager_ == true) {
                    TlMemManager& rMemManager = TlMemManager::getInstance();
                    rMemManager.deallocate((char*)this->data_[i]);
                } else {
                    delete[] this->data_[i];
                }
                this->data_[i] = NULL;
            }

            this->data_[i] = pNew;
        }
    }
}


void TlRowVectorMatrix2::set(index_type row, index_type col, double value)
{
    const div_t turns = std::div(row, this->allProcs_);
    if (turns.rem == this->rank_) {
        const index_type row = turns.quot;
        assert(row < this->numOfLocalRows_);
        assert(col < this->getNumOfCols());
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
        assert(this->getNumOfCols() <= this->reserveCols_);
        answer = TlVector(this->data_[row], this->getNumOfCols());
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
        std::copy(this->data_[row],
                  this->data_[row] + copySize,
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

    const index_type numOfLocalRows = this->numOfLocalRows_;
    for (index_type r = 0; r < numOfLocalRows; ++r) {
        const index_type row = r * this->allProcs_ + this->rank_;
        for (index_type col = 0; col < numOfCols; ++col) {
            answer.set(row, col, this->data_[r][col]);
        }
    }
    
    return answer;
}


