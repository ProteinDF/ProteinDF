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
#include "TlListMatrix.h"

TlListMatrix::TlListMatrix(const std::size_t reserveSize) : size_(0)
{
    this->elements_.clear();
}

TlListMatrix::~TlListMatrix()
{
}

void TlListMatrix::clear()
{
    this->size_ = 0;
    this->elements_.clear();
}

void TlListMatrix::add(const int row, const int col, const double value)
{
    this->elements_.push_back(Element(row, col, value));
    ++this->size_;
}

TlListMatrix::Element TlListMatrix::pop()
{
    Element answer(-1, -1, 0.0);
    if (this->size() > 0) {
        answer = this->elements_.front();
        this->elements_.pop_front();
        --this->size_;
    }

    return answer;
}
