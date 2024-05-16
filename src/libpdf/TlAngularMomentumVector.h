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

#ifndef TLANGULARMOMENTUMVECTOR_H
#define TLANGULARMOMENTUMVECTOR_H

#include <cassert>
#include <iostream>

#include "TlUtils.h"

/// @class angular momentum vector
///   a = (a_x, a_y, a_z)
///   同じAngular Momentumのvector数は(L+1)*(L+2)/2
class TlAngularMomentumVector {
private:
    enum { X = 0,
           Y = 1,
           Z = 2 };

public:
    TlAngularMomentumVector(int x = 0, int y = 0, int z = 0);

    TlAngularMomentumVector(const TlAngularMomentumVector& rhs);

    /// angular momentumを返す
    int angularMomentum() const;

    int index() const;

    bool isExist() const;

    /// 指定された軸[x(=0), y(=1), z(=3)] の要素を返す
    int get(const int i) const;

    TlAngularMomentumVector operator+=(const TlAngularMomentumVector& rhs);

    TlAngularMomentumVector operator-=(const TlAngularMomentumVector& rhs);

    friend TlAngularMomentumVector operator+(
        const TlAngularMomentumVector& rhs1,
        const TlAngularMomentumVector& rhs2);

    friend TlAngularMomentumVector operator-(
        const TlAngularMomentumVector& rhs1,
        const TlAngularMomentumVector& rhs2);

    std::string debugOut() const;

private:
    int v_[3];
};

inline bool TlAngularMomentumVector::isExist() const {
    // return ((v_[0] >= 0) && (v_[1] >= 0) && (v_[2] >= 0));
    return ((this->v_[0] | this->v_[1]) | this->v_[2]) >= 0;
}

inline int TlAngularMomentumVector::index() const {
    // (y * (y + 1) + z * (2*y + z + 3)) /2 = ((y+z)*(y+z+1)+2z)/2

    // const int y = this->v_[Y];
    const int z = this->v_[Z];
    // const int yz = this->v_[Y] + z;
    const int yz = this->v_[Y] + this->v_[Z];
    return ((yz * (yz + 1)) >> 1) + z;
}

#endif  // TLANGULARMOMENTUMVECTOR_H
