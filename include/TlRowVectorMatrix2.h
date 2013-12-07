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

#ifndef TLROWVECTORMATRIX2_H
#define TLROWVECTORMATRIX2_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include "TlMatrixObject.h"
#include "TlMatrix.h"

class TlRowVectorMatrix2 {
public:
    typedef TlMatrixObject::index_type index_type;

public:
    explicit TlRowVectorMatrix2(index_type row = 1, index_type col = 1,
                                int allProcs = 1, int rank = 0,
                                bool isUsingMemManager = false);
    TlRowVectorMatrix2(const TlRowVectorMatrix2& rhs);

    ~TlRowVectorMatrix2();
        
public:

    void resize(index_type row, index_type col);
    index_type getNumOfRows() const {
        return this->numOfRows_;
    };

    index_type getNumOfCols() const {
        return this->numOfCols_;
    };

    index_type getNumOfLocalRows() const {
        return this->numOfLocalRows_;
    }

    /// 列数のcapacityを設定する
    void reserve_cols(index_type col);

    void set(index_type row, index_type col, double value);
        
    TlVector getRowVector(index_type row) const;
    index_type getRowVector(index_type row, double *pBuf, index_type maxColSize) const;
        
    int getPEinChargeByRow(index_type row) const;
        
    TlMatrix getTlMatrix() const;

private:
    TlRowVectorMatrix2& operator=(const TlRowVectorMatrix2& rhs);

private:
    index_type numOfRows_;
    index_type numOfCols_;
    index_type reserveCols_; // 列数のメモリ確保量

    int allProcs_;
    int rank_;

    index_type numOfLocalRows_;
    std::vector<double* > data_;

    bool isUsingMemManager_;
};

#endif // TLROWVECTORMATRIX2_H
