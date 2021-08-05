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

#include "tl_matrix_object.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <string>

#include "TlUtils.h"

std::string TlMatrixObject::matrixTypeStr(TlMatrixObject::MatrixType matrixType) {
    static std::map<TlMatrixObject::MatrixType, std::string> matrixType2str;
    if (matrixType2str.size() == 0) {
        // initialize
        matrixType2str[UNDEFINED] = "UNDEFINED";
        matrixType2str[RSFD] = "RSFD";
        matrixType2str[CSFD] = "CSFD";
        matrixType2str[RLHD] = "RLHD";
        matrixType2str[RUHD] = "RUHD";
        matrixType2str[CLHD] = "CLHD";
        matrixType2str[CUHD] = "CUHD";
        matrixType2str[COOF] = "COOF";
        matrixType2str[COOS] = "COOS";
        matrixType2str[ABGD] = "ABGD";
    }

    std::map<TlMatrixObject::MatrixType, std::string>::const_iterator it = matrixType2str.find(matrixType);
    std::string answer = "";
    if (it != matrixType2str.end()) {
        answer = it->second;
    }

    return answer;
}

TlMatrixObject::TlMatrixObject(const MatrixType matrixType, const index_type row, const index_type col)
    : matrixType_(matrixType), row_(row), col_(col), log_(TlLogging::getInstance()) {
}

TlMatrixObject::size_type TlMatrixObject::getNumOfElements() const {
    TlMatrixObject::size_type elements = 0;

    switch (this->matrixType_) {
        case TlMatrixObject::CSFD:
            elements = this->getNumOfElements_CSFD();
            break;

        case TlMatrixObject::RSFD:
            elements = this->getNumOfElements_RSFD();
            break;

        case TlMatrixObject::RLHD:
            elements = this->getNumOfElements_RLHD();
            break;

        default:
            this->log_.critical(TlUtils::format("program error: unknown matrix type @%s, %d", __FILE__, __LINE__));
            break;
    }

    return elements;
}

TlMatrixObject::size_type TlMatrixObject::getNumOfElements_CSFD() const {
    const size_type row = static_cast<size_type>(this->getNumOfRows());
    const size_type col = static_cast<size_type>(this->getNumOfCols());

    const size_type elements = row * col;
    return elements;
}

TlMatrixObject::size_type TlMatrixObject::getNumOfElements_RSFD() const {
    // the same value of CSFD!
    return this->getNumOfElements_CSFD();
}

TlMatrixObject::size_type TlMatrixObject::getNumOfElements_RLHD() const {
    const size_type dim = static_cast<size_type>(this->getNumOfRows());
    assert(dim == static_cast<size_type>(this->getNumOfCols()));

    const size_type elements = dim * (dim + 1) / 2;
    return elements;
}

TlMatrixObject::size_type TlMatrixObject::getIndex(index_type row, index_type col) const {
    TlMatrixObject::size_type index = 0;

    switch (this->matrixType_) {
        case TlMatrixObject::CSFD:
            index = this->getIndex_CSFD(row, col);
            break;

        case TlMatrixObject::RSFD:
            index = this->getIndex_RSFD(row, col);
            break;

        case TlMatrixObject::RLHD:
            index = this->getIndex_RLHD(row, col);
            break;

        default:
            this->log_.critical(TlUtils::format("program error: unknown matrix type @%s, %d", __FILE__, __LINE__));
            break;
    }

    return index;
}

TlMatrixObject::size_type TlMatrixObject::getIndex_RSFD(index_type row, index_type col) const {
    assert(0 <= row);
    assert(row < this->getNumOfRows());
    assert(0 <= col);
    assert(col < this->getNumOfCols());

    const size_type r = static_cast<size_type>(row);
    const size_type c = static_cast<size_type>(col);
    const size_type maxC = static_cast<size_type>(this->getNumOfCols());

    const size_type addr = r * maxC + c;
    return addr;
}

TlMatrixObject::size_type TlMatrixObject::getIndex_CSFD(index_type row, index_type col) const {
    assert(0 <= row);
    assert(row < this->getNumOfRows());
    assert(0 <= col);
    assert(col < this->getNumOfCols());

    const size_type r = static_cast<size_type>(row);
    const size_type c = static_cast<size_type>(col);
    const size_type maxR = static_cast<size_type>(this->getNumOfRows());

    const size_type addr = c * maxR + r;
    return addr;
}

TlMatrixObject::size_type TlMatrixObject::getIndex_RLHD(index_type row, index_type col) const {
    if (row < col) {
        std::swap(row, col);
    }

    const size_type r = static_cast<size_type>(row);
    const size_type c = static_cast<size_type>(col);

    const size_type addr = r * (r + 1) / 2 + c;
    return addr;
}

// TlVector TlMatrixObject::getRowVector(const int row) const {
//   assert((0 <= row) && (row < this->getNumOfRows()));
//
//   const int numOfCols = this->getNumOfCols();
//   TlVector answer(numOfCols);
//
//   for (int i = 0; i < numOfCols; ++i) {
//     answer[i] = this->get(row, i);
//   }
//
//   return answer;
// }
//
// TlVector TlMatrixObject::getColVector(const int col) const {
//   assert((0 <= col) && (col < this->getNumOfCols()));
//
//   const int numOfRows = this->getNumOfRows();
//   TlVector answer(numOfRows);
//
//   for (int i = 0; i < numOfRows; ++i) {
//     answer[i] = this->get(i, col);
//   }
//
//   return answer;
// }

std::ostream& operator<<(std::ostream& stream, const TlMatrixObject& mat) {
    const TlMatrixObject::index_type numOfRows = mat.getNumOfRows();
    const TlMatrixObject::index_type numOfCols = mat.getNumOfCols();

    for (TlMatrixObject::index_type ord = 0; ord < numOfCols; ord += 10) {
        stream << "       ";
        for (TlMatrixObject::index_type j = ord; ((j < ord + 10) && (j < numOfCols)); ++j) {
            stream << TlUtils::format("   %5d th", j + 1);
        }
        stream << "\n ----";

        for (TlMatrixObject::index_type j = ord; ((j < ord + 10) && (j < numOfCols)); ++j) {
            stream << "-----------";
        }
        stream << "----\n";

        for (TlMatrixObject::index_type i = 0; i < numOfRows; ++i) {
            stream << TlUtils::format(" %5d  ", i + 1);

            for (TlMatrixObject::index_type j = ord; ((j < ord + 10) && (j < numOfCols)); ++j) {
                stream << TlUtils::format(" %10.6lf", mat.get(i, j));
            }
            stream << "\n";
        }
        stream << "\n\n";
    }

    return stream;
}
