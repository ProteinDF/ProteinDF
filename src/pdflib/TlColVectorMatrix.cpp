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
#include "TlColVectorMatrix.h"

TlColVectorMatrix::TlColVectorMatrix(const index_type row,
                                     const index_type col,
                                     const int numOfSubunits,
                                     const int subunitID,
                                     bool isUsingMemManager)
    : TlVectorMatrixObject(col, row, numOfSubunits, subunitID, isUsingMemManager) {
}


TlColVectorMatrix::~TlColVectorMatrix()
{
}


void TlColVectorMatrix::resize(const index_type row,
                               const index_type col)
{
    TlVectorMatrixObject::resize(col, row);
}


void TlColVectorMatrix::reserveRowSize(const index_type reserveRowSize)
{
    TlVectorMatrixObject::reserveVectorSize(reserveRowSize);
}

void TlColVectorMatrix::set(const index_type row,
                            const index_type col, 
                            const double value)
{
    TlVectorMatrixObject::set(col, row, value);
}


double TlColVectorMatrix::get(const index_type row,
                              const index_type col) const
{
    return TlVectorMatrixObject::get(col, row);
}


TlMatrix TlColVectorMatrix::getTlMatrixObject() const
{
    const index_type numOfRows = this->getNumOfRows();
    const index_type numOfCols = this->getNumOfCols();
    TlMatrix answer(numOfRows, numOfCols);

    std::vector<double> v;
    for (index_type c = 0; c < numOfCols; ++c) {
        v = TlVectorMatrixObject::getVector(c);
        assert(v.size() == numOfRows);
        for (index_type r = 0; r < numOfRows; ++r) {
            answer.set(r, c, v[r]);
        }
    }
    
    return answer;
}

