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

TlColVectorMatrix::TlColVectorMatrix(const TlMatrix& rhs,
                                     const int numOfSubunits,
                                     const int subunitID,
                                     bool isUsingMemManager)
    : TlVectorMatrixObject(rhs.getNumOfCols(), rhs.getNumOfRows(),
                           numOfSubunits, subunitID, isUsingMemManager)
{
    const index_type numOfRows = rhs.getNumOfRows();
    const index_type numOfCols = rhs.getNumOfCols();
    for (index_type r = 0; r < numOfRows; ++r) {
        for (index_type c = 0; c < numOfCols; ++c) {
            this->set(r, c, rhs.get(r, c));
        }
    }
}

TlColVectorMatrix::TlColVectorMatrix(const TlColVectorMatrix& rhs)
    : TlVectorMatrixObject(rhs)
{
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
    TlVectorMatrixObject::set_to_vm(col, row, value);
}


void TlColVectorMatrix::add(const index_type row,
                            const index_type col, 
                            const double value)
{
    TlVectorMatrixObject::add_to_vm(col, row, value);
}


double TlColVectorMatrix::get(const index_type row,
                              const index_type col) const
{
    return TlVectorMatrixObject::get_from_vm(col, row);
}

TlVector TlColVectorMatrix::getColVector(const index_type col) const
{
    return TlVectorMatrixObject::getVector(col);
}

TlMatrix TlColVectorMatrix::getTlMatrixObject() const
{
    const index_type numOfRows = this->getNumOfRows();
    const index_type numOfCols = this->getNumOfCols();
    TlMatrix answer(numOfRows, numOfCols);

    for (index_type c = 0; c < numOfCols; ++c) {
        const TlVector v = TlVectorMatrixObject::getVector(c);
        assert(v.getSize() == numOfRows);
        for (index_type r = 0; r < numOfRows; ++r) {
            answer.set(r, c, v[r]);
        }
    }
    
    return answer;
}

void TlColVectorMatrix::saveByTlRowVectorMatrix(const std::string& basename) const
{
    TlVectorMatrixObject::saveByTheOtherType(basename);
}


