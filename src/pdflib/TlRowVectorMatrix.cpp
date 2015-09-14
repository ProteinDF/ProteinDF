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
#include "TlRowVectorMatrix.h"

TlRowVectorMatrix::TlRowVectorMatrix(const index_type row,
                                     const index_type col,
                                     const int numOfSubunits,
                                     const int subunitID,
                                     bool isUsingMemManager)
    : TlVectorMatrixObject(row, col, numOfSubunits, subunitID, isUsingMemManager) {
}

TlRowVectorMatrix::TlRowVectorMatrix(const TlMatrix& rhs, 
                                     const int numOfSubunits,
                                     const int subunitID,
                                     bool isUsingMemManager)
    : TlVectorMatrixObject(rhs.getNumOfRows(), rhs.getNumOfCols(),
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

TlRowVectorMatrix::TlRowVectorMatrix(const TlRowVectorMatrix& rhs)
    : TlVectorMatrixObject(rhs)
{
}

TlRowVectorMatrix::~TlRowVectorMatrix()
{
}


void TlRowVectorMatrix::resize(const index_type row,
                               const index_type col)
{
    TlVectorMatrixObject::resize(row, col);
}


void TlRowVectorMatrix::reserveColSize(const index_type reserveColSize)
{
    TlVectorMatrixObject::reserveVectorSize(reserveColSize);
}


void TlRowVectorMatrix::set(const index_type row,
                            const index_type col, 
                            const double value)
{
    TlVectorMatrixObject::set_to_vm(row, col, value);
}

void TlRowVectorMatrix::add(const index_type row,
                            const index_type col, 
                            const double value)
{
    TlVectorMatrixObject::add_to_vm(row, col, value);
}


double TlRowVectorMatrix::get(const index_type row,
                              const index_type col) const
{
    return TlVectorMatrixObject::get_from_vm(row, col);
}

TlVector TlRowVectorMatrix::getRowVector(const index_type row) const
{
    return TlVectorMatrixObject::getVector(row);
}

TlMatrix TlRowVectorMatrix::getTlMatrixObject() const
{
    const index_type numOfRows = this->getNumOfRows();
    const index_type numOfCols = this->getNumOfCols();
    TlMatrix answer(numOfRows, numOfCols);

    for (index_type r = 0; r < numOfRows; ++r) {
        TlVector v = TlVectorMatrixObject::getVector(r);
        assert(v.getSize() == numOfCols);
        for (index_type c = 0; c < numOfCols; ++c) {
            answer.set(r, c, v[c]);
        }
    }
    
    return answer;
}

void TlRowVectorMatrix::saveByTlColVectorMatrix(const std::string& basename) const
{
    TlVectorMatrixObject::saveByTheOtherType(basename);
}
