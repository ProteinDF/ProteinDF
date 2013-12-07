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

#ifndef TLANGULARMOMENTUMVECTORSET_H
#define TLANGULARMOMENTUMVECTORSET_H

#include <vector>
#include <bitset>
#include "TlAngularMomentumVector.h"

#define MAX_ANGULAR_MOMENTUM (17) // 0~16までOK

/// Angular Momentum Vector Set
///
/// angular momentum
class TlAngularMomentumVectorSet
{
public:
    TlAngularMomentumVectorSet(const int angularMomentum);

    /// angular momentum vector の総数を返す
    int size() const;

    /// 指定されたインデックスのangular momentum vectorを返す
    /// @ref count
    const TlAngularMomentumVector& get(const int index) const;

private:
    struct AMV_DB {
    public:
        AMV_DB() {
            this->data.resize(MAX_ANGULAR_MOMENTUM);
            for (int amv = 0; amv < MAX_ANGULAR_MOMENTUM; ++amv) {
                std::vector<TlAngularMomentumVector> temp;
                const int numOfAMVs = (amv + 1) * (amv + 2) / 2;
                temp.resize(numOfAMVs);
                
                int index = 0;
                for (char x = amv; x >= 0; --x) {
                    for (char y = (amv - x); y >= 0; --y) {
                        const char z = amv - (x + y);
                        temp[index] = TlAngularMomentumVector(x, y, z);
                        ++index;
                    }
                }
                this->data[amv] = temp;
            }
        };

    public:
        std::vector<std::vector<TlAngularMomentumVector> > data;
    };
    
private:
    static const AMV_DB db_;
    int angularMomentum_;
};


inline TlAngularMomentumVectorSet::TlAngularMomentumVectorSet(const int angularMomentum)
    : angularMomentum_(angularMomentum) {
    assert((0 <= angularMomentum) && (angularMomentum < MAX_ANGULAR_MOMENTUM));
}

inline int TlAngularMomentumVectorSet::size() const
{
    return (angularMomentum_+1)*(angularMomentum_+2) / 2;
}

inline const TlAngularMomentumVector&
TlAngularMomentumVectorSet::get(const int index) const
{
    assert((0 <= index) && (index < this->size()));
    return this->db_.data[this->angularMomentum_][index];
}


#endif // TLANGULARMOMENTUMVECTORSET_H
