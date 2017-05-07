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

#ifndef TLPARTIALSYMMETRICMATRIX_H
#define TLPARTIALSYMMETRICMATRIX_H

#include <vector>
#include "TlUtils.h"
#include "TlPartialMatrix.h"

class TlPartialSymmetricMatrix : public TlPartialMatrix {
public:
    explicit TlPartialSymmetricMatrix(index_type globalSize =0,
                                      index_type startRow =0,
                                      index_type startCol =0,
                                      index_type range =0);
    TlPartialSymmetricMatrix(index_type globalSize,
                             index_type startRow,
                             index_type startCol,
                             index_type rowRange,
                             index_type colRange);
    TlPartialSymmetricMatrix(const TlPartialSymmetricMatrix& rhs);
    virtual ~TlPartialSymmetricMatrix();

public:
    TlPartialSymmetricMatrix& operator=(const TlPartialSymmetricMatrix& rhs);

public:
    virtual void set(index_type globalRow, index_type globalCol, double value);
    virtual void add(index_type globalRow, index_type globalCol, double value);
    virtual double get(int globalRow, index_type globalCol) const;
    virtual double getLocal(index_type localRow, index_type localCol) const;

    virtual TlVector getRowVector(index_type row) const;
    virtual TlVector getColVector(index_type col) const;

public:
    template <typename T>
    void print(T& out) const;

protected:
    virtual size_type index(index_type globalRow, index_type globalCol) const;

private:
    //
    friend class TlCommunicate;
};


// =============================================================================
inline TlMatrixObject::size_type TlPartialSymmetricMatrix::index(index_type globalRow, index_type globalCol) const
{
    // 行・列の優先順位をつけるために必要
    if (globalRow < globalCol) {
        std::swap(globalRow, globalCol);
    }
    int answer = TlPartialMatrix::index(globalRow, globalCol);
    if (answer == -1) {
        answer = TlPartialMatrix::index(globalCol, globalRow);
    }

    return answer;
}


template <typename T>
void TlPartialSymmetricMatrix::print(T& out) const
{
    out << "start row = " << this->startRow_ << "\n";
    out << "start col = " << this->startCol_ << "\n";
    const int nNumOfRows = this->localNumOfRows_;
    const int nNumOfCols = this->localNumOfCols_;
//   if (this->startRow_ == this->startCol_) {
//     const int nNumOfDim = nNumOfRows;
//     out << "\n\n";
//     for (int ord = 0; ord < nNumOfDim; ord += 10) {
//       out << "       ";
//       for (int j = ord; ((j < ord+10) && (j < nNumOfDim)); ++j) {
//  out << TlUtils::format("   %5d th", j+1);
//       }
//       out << "\n" << " ----";

//       for (int j = ord; ((j < ord+10) && (j < nNumOfDim)); ++j) {
//  out << "-----------";
//       }
//       out << "----\n";

//       for (int i = 0; i < nNumOfDim; ++i) {
//  out << TlUtils::format(" %5d  ", i+1);

//  for (int j = ord; ((j < ord+10) && (j < nNumOfDim)); ++j) {
//    if (j > i){
//      out << "    ----   ";
//    } else {
//      out << TlUtils::format(" %10.6lf", (*this)(i,j));
//    }
//  }
//  out << "\n";
//       }
//       out << "\n\n";
//     }
//     out.flush();
//   } else {
    for (int ord = 0; ord < nNumOfCols; ord += 10) {
        out << "       ";
        for (int j = ord; ((j < ord+10) && (j < nNumOfCols)); ++j) {
            out << TlUtils::format("   %5d th", j+1);
        }
        out << "\n ----";

        for (int j = ord; ((j < ord+10) && (j < nNumOfCols)); ++j) {
            out << "-----------";
        }
        out << "----\n";

        for (int i = 0; i < nNumOfRows; ++i) {
            out << TlUtils::format(" %5d  ", i+1);

            for (int j = ord; ((j < ord+10) && (j < nNumOfCols)); ++j) {
                out << TlUtils::format(" %10.6lf", this->getLocal(i, j));
            }
            out << "\n";
        }
        out << "\n\n";
    }
    out.flush();
//   }
}

#endif // TLPARTIALSYMMETRICMATRIX_H

