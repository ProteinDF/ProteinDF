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

#ifndef TLLISTMATRIX_H
#define TLLISTMATRIX_H

#include <deque>

class TlListMatrix {
public:
    struct Element {
public:
        Element(int r = 0, int c = 0, double v = 0.0) : row(r), col(c), value(v) {
        }

public:
        int row;
        int col;
        double value;
    };

public:
    TlListMatrix(std::size_t reserveSize = 0);
    ~TlListMatrix();

public:
    void clear();

    std::size_t size() const {
        return this->size_;
    };

    void add(int row, int col, double value);

    Element pop();

private:
    std::size_t size_;
    std::deque<Element> elements_;
};

#endif // TLLISTMATRIX_H
