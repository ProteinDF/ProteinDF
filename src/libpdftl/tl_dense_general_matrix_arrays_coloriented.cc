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

#include <cassert>
#include <iostream>

#include "TlUtils.h"
#include "tl_dense_general_matrix_arrays_coloriented.h"

TlDenseGeneralMatrix_arrays_ColOriented::
    TlDenseGeneralMatrix_arrays_ColOriented(const index_type row,
                                            const index_type col,
                                            const int numOfSubunits,
                                            const int subunitID,
                                            bool isUsingMemManager)
    : TlDenseMatrix_arrays_Object(col, row, numOfSubunits, subunitID,
                                  isUsingMemManager) {}

TlDenseGeneralMatrix_arrays_ColOriented::
    TlDenseGeneralMatrix_arrays_ColOriented(
        const TlDenseGeneralMatrix_Lapack& rhs, const int numOfSubunits,
        const int subunitID, bool isUsingMemManager)
    : TlDenseMatrix_arrays_Object(rhs.getNumOfCols(), rhs.getNumOfRows(),
                                  numOfSubunits, subunitID, isUsingMemManager) {
    const index_type numOfRows = rhs.getNumOfRows();
    const index_type numOfCols = rhs.getNumOfCols();
    for (index_type r = 0; r < numOfRows; ++r) {
        for (index_type c = 0; c < numOfCols; ++c) {
            this->set(r, c, rhs.get(r, c));
        }
    }
}

TlDenseGeneralMatrix_arrays_ColOriented::
    TlDenseGeneralMatrix_arrays_ColOriented(
        const TlDenseGeneralMatrix_arrays_ColOriented& rhs)
    : TlDenseMatrix_arrays_Object(rhs) {}

TlDenseGeneralMatrix_arrays_ColOriented::
    TlDenseGeneralMatrix_arrays_ColOriented(const std::string& path)
    : TlDenseMatrix_arrays_Object() {
    this->load(path);
}

TlDenseGeneralMatrix_arrays_ColOriented::
    ~TlDenseGeneralMatrix_arrays_ColOriented() {}

void TlDenseGeneralMatrix_arrays_ColOriented::resize(const index_type row,
                                                     const index_type col) {
    TlDenseMatrix_arrays_Object::resize(col, row);
}

void TlDenseGeneralMatrix_arrays_ColOriented::reserveRowSize(
    const index_type reserveRowSize) {
    TlDenseMatrix_arrays_Object::reserveVectorSize(reserveRowSize);
}

void TlDenseGeneralMatrix_arrays_ColOriented::set(const index_type row,
                                                  const index_type col,
                                                  const double value) {
    TlDenseMatrix_arrays_Object::set_to_vm(col, row, value);
}

void TlDenseGeneralMatrix_arrays_ColOriented::add(const index_type row,
                                                  const index_type col,
                                                  const double value) {
    TlDenseMatrix_arrays_Object::add_to_vm(col, row, value);
}

double TlDenseGeneralMatrix_arrays_ColOriented::get(
    const index_type row, const index_type col) const {
    return TlDenseMatrix_arrays_Object::get_from_vm(col, row);
}

std::vector<double> TlDenseGeneralMatrix_arrays_ColOriented::getColVector(
    const index_type col) const {
    return TlDenseMatrix_arrays_Object::getVector(col);
}

void TlDenseGeneralMatrix_arrays_ColOriented::getColVector(
    const index_type col, double* pBuf, const index_type length) const {
    TlDenseMatrix_arrays_Object::getVector(col, pBuf, length);
}

TlDenseGeneralMatrix_Lapack
TlDenseGeneralMatrix_arrays_ColOriented::getTlMatrixObject() const {
    const index_type numOfRows = this->getNumOfRows();
    const index_type numOfCols = this->getNumOfCols();
    TlDenseGeneralMatrix_Lapack answer(numOfRows, numOfCols);

    for (index_type c = 0; c < numOfCols; ++c) {
        const TlDenseVector_Lapack v =
            TlDenseMatrix_arrays_Object::getVector(c);
        assert(v.getSize() == numOfRows);
        for (index_type r = 0; r < numOfRows; ++r) {
            answer.set(r, c, v.get(r));
        }
    }

    return answer;
}

// -----------------------------------------------------------------------------
std::ostream& operator<<(std::ostream& stream,
                         const TlDenseGeneralMatrix_arrays_ColOriented& mat) {
    const TlMatrixObject::index_type numOfRows = mat.getNumOfRows();
    const TlMatrixObject::index_type numOfCols = mat.getNumOfCols();

    for (TlMatrixObject::index_type ord = 0; ord < numOfCols; ord += 10) {
        stream << "       ";
        for (TlMatrixObject::index_type j = ord;
             ((j < ord + 10) && (j < numOfCols)); ++j) {
            stream << TlUtils::format("   %5d th", j + 1);
        }
        stream << "\n ----";

        for (TlMatrixObject::index_type j = ord;
             ((j < ord + 10) && (j < numOfCols)); ++j) {
            stream << "-----------";
        }
        stream << "----\n";

        for (TlMatrixObject::index_type i = 0; i < numOfRows; ++i) {
            stream << TlUtils::format(" %5d  ", i + 1);

            for (TlMatrixObject::index_type j = ord;
                 ((j < ord + 10) && (j < numOfCols)); ++j) {
                stream << TlUtils::format(" %10.6lf", mat.get(i, j));
            }
            stream << "\n";
        }
        stream << "\n\n";
    }

    return stream;
}

// -----------------------------------------------------------------------------
void TlDenseGeneralMatrix_arrays_ColOriented::
    saveByTlDenseGeneralMatrix_arrays_RowOriented(
        const std::string& basename) const {
    TlDenseMatrix_arrays_Object::saveByTheOtherType(basename);
}
