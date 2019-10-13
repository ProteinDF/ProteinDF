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

#include "TlAngularMomentumVector.h"
#include "TlUtils.h"

TlAngularMomentumVector::TlAngularMomentumVector(char x, char y, char z) {
    this->v_[X] = x;
    this->v_[Y] = y;
    this->v_[Z] = z;
}

TlAngularMomentumVector::TlAngularMomentumVector(
    const TlAngularMomentumVector& rhs) {
    this->v_[X] = rhs.v_[X];
    this->v_[Y] = rhs.v_[Y];
    this->v_[Z] = rhs.v_[Z];
}

int TlAngularMomentumVector::angularMomentum() const {
    return this->v_[X] + this->v_[Y] + this->v_[Z];
}

int TlAngularMomentumVector::get(const int i) const {
    assert((0 <= i) && (i < 3));
    return this->v_[i];
}

TlAngularMomentumVector TlAngularMomentumVector::operator+=(
    const TlAngularMomentumVector& rhs) {
    this->v_[X] += rhs.v_[X];
    this->v_[Y] += rhs.v_[Y];
    this->v_[Z] += rhs.v_[Z];

    return *this;
}

TlAngularMomentumVector TlAngularMomentumVector::operator-=(
    const TlAngularMomentumVector& rhs) {
    this->v_[X] -= rhs.v_[X];
    this->v_[Y] -= rhs.v_[Y];
    this->v_[Z] -= rhs.v_[Z];

    return *this;
}

TlAngularMomentumVector operator+(const TlAngularMomentumVector& rhs1,
                                  const TlAngularMomentumVector& rhs2) {
    TlAngularMomentumVector ans = rhs1;
    ans += rhs2;
    return ans;
}

TlAngularMomentumVector operator-(const TlAngularMomentumVector& rhs1,
                                  const TlAngularMomentumVector& rhs2) {
    TlAngularMomentumVector ans = rhs1;
    ans -= rhs2;
    return ans;
}

std::string TlAngularMomentumVector::debugOut() const {
    return TlUtils::format("[%d %d %d]", v_[0], v_[1], v_[2]);
}
