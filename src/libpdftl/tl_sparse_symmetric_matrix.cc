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

#include "tl_sparse_symmetric_matrix.h"

TlSparseSymmetricMatrix::TlSparseSymmetricMatrix(const int size)
    : TlSparseMatrix(size, size) {}

TlSparseSymmetricMatrix::TlSparseSymmetricMatrix(
    const TlSparseSymmetricMatrix& rhs)
    : TlSparseMatrix(rhs) {}

TlSparseSymmetricMatrix::~TlSparseSymmetricMatrix() {}

void TlSparseSymmetricMatrix::resize(const int size) {
    TlSparseMatrix::resize(size, size);
}

TlDenseVector_Lapack TlSparseSymmetricMatrix::getRowVector(
    const index_type row) const {
    const int numOfCols = this->getNumOfCols();
    TlDenseVector_Lapack ans(numOfCols);

    if (this->m_aMatrix.size() > (std::size_t)(numOfCols / 2)) {
        for (index_type i = 0; i < numOfCols; ++i) {
            ans.set(i, this->get(row, i));
        }
    } else {
        const_iterator pEnd = this->end();
        for (const_iterator p = this->begin(); p != pEnd; ++p) {
            const int r = p->first.row;
            const int c = p->first.col;
            if (r == row) {
                ans.set(c, p->second);
            } else if (c == row) {
                ans.set(r, p->second);
            }
        }
    }

    return ans;
}

TlDenseVector_Lapack TlSparseSymmetricMatrix::getColVector(
    const index_type col) const {
    const int numOfRows = this->getNumOfRows();
    TlDenseVector_Lapack ans(numOfRows);

    if (this->m_aMatrix.size() > (std::size_t)(numOfRows / 2)) {
        for (index_type i = 0; i < numOfRows; ++i) {
            ans.set(i, this->get(i, col));
        }
    } else {
        const_iterator pEnd = this->end();
        for (const_iterator p = this->begin(); p != pEnd; ++p) {
            const int r = p->first.row;
            const int c = p->first.col;
            if (c == col) {
                ans.set(r, p->second);
            } else if (r == col) {
                ans.set(c, p->second);
            }
        }
    }

    return ans;
}

// const TlSparseSymmetricMatrix& TlSparseSymmetricMatrix::dot(const
// TlSparseSymmetricMatrix& X)
// {
//     iterator p = this->begin();
//     iterator pEnd = this->end();
//     const_iterator q = X.begin();
//     const_iterator qEnd = X.end();

//     while ((p != pEnd) && (q != qEnd)) {
//         const unsigned long p_index = p->first;
//         const unsigned long q_index = q->first;
//         if (p_index < q_index) {
//             ++p;
//         } else if (p_index > q_index) {
//             ++q;
//         } else {
//             // p_index == q_index
//             this->m_aMatrix[p_index] *= X.m_aMatrix[p_index];
//             ++p;
//             ++q;
//         }
//     }

//     return (*this);
// }
