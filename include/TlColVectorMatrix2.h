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

#ifndef TLCOLVECTORMATRIX2_H
#define TLCOLVECTORMATRIX2_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include "TlMatrixObject.h"
#include "TlMatrix.h"
#include "TlLogging.h"

class TlColVectorMatrix2 {
public:
    typedef TlMatrixObject::index_type index_type;

public:
    explicit TlColVectorMatrix2(index_type row = 1, index_type col = 1,
                                int allProcs = 1, int rank = 0,
                                bool isUsingMemManager = false);
    ~TlColVectorMatrix2();
        
public:

    void resize(index_type row, index_type col);
    index_type getNumOfRows() const {
        return this->numOfRows_;
    };

    index_type getNumOfCols() const {
        return this->numOfCols_;
    };

    /// 行数のcapacityを設定する
    void reserve_rows(index_type row);

    void set(index_type row, index_type col, double value);
        
    TlVector getColVector(index_type col) const;
    index_type getColVector(index_type col, double *pBuf, index_type maxRowSize) const;
        
    int getPEinChargeByCol(index_type col) const;

    TlMatrix getTlMatrix() const;
        
public:
    int getNumOfAllProcs() const {
        return this->allProcs_;
    };

    int getRank() const {
        return this->rank_;
    };

    void save(const std::string& basename) const;
    void load(const std::string& basename);

private:
    TlColVectorMatrix2(const TlColVectorMatrix2& rhs);
    TlColVectorMatrix2& operator=(const TlColVectorMatrix2& rhs);

private:
    TlLogging& log_;

    index_type numOfRows_;
    index_type numOfCols_;
    index_type reserveRows_; // 行数のメモリ確保量

    int allProcs_;
    int rank_;

    index_type numOfLocalCols_;
    std::vector<double* > data_;

    bool isUsingMemManager_;
};

#endif // TLCOLVECTORMATRIX2_H
