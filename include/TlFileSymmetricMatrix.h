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

#ifndef TLFILESYMMETRICMATRIX_H
#define TLFILESYMMETRICMATRIX_H

#include "TlFileMatrix.h"
#include "TlPartialSymmetricMatrix.h"

class TlFileSymmetricMatrix : public TlFileMatrix {
public:
    explicit TlFileSymmetricMatrix(const std::string& filePath, int dim = 0, size_t cacheSize = DEFAULT_CACHE_SIZE);
    virtual ~TlFileSymmetricMatrix();

public:
    void add(int row, int col, double value);
    void add(const TlPartialSymmetricMatrix& psm);

    virtual TlFileSymmetricMatrix& operator*=(double coef);
    virtual TlFileSymmetricMatrix& operator/=(double coef) {
        return this->operator*=(1.0 / coef);
    }

    TlPartialSymmetricMatrix getPartialMatrix(const int startRow, const int startCol, const int range) const;

protected:
    virtual void open();

protected:
    virtual bool readHeader();

    virtual size_t index(int row, int col) const {
        if (row < col) {
            std::swap(row, col);
        }

        return (std::size_t(col) + (std::size_t(row +1) * std::size_t(row)) / std::size_t(2));
    }

    virtual size_t maxIndex() const {
        return (std::size_t(this->getNumOfRows()) * std::size_t(this->getNumOfCols() +1) / std::size_t(2));
    }
};

#endif // TLFILESYMMETRICMATRIX_H
