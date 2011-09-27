#ifndef TLANGULARMOMENTUMVECTOR_H
#define TLANGULARMOMENTUMVECTOR_H

#include <iostream>
#include <cassert>
#include "TlUtils.h"

/// @class angular momentum vector
///   a = (a_x, a_y, a_z)
///   同じAngular Momentumのvector数は(L+1)*(L+2)/2
class TlAngularMomentumVector
{
private:
    enum {
        X = 0,
        Y = 1,
        Z = 2
    };
    
public:
    TlAngularMomentumVector(char x = 0,
                            char y = 0,
                            char z = 0) {
        this->v_[X] = x;
        this->v_[Y] = y;
        this->v_[Z] = z;
    }

    TlAngularMomentumVector(const TlAngularMomentumVector& rhs) {
        this->v_[X] = rhs.v_[X];
        this->v_[Y] = rhs.v_[Y];
        this->v_[Z] = rhs.v_[Z];
    }

    /// angular momentumを返す
    unsigned int angularMomentum() const {
        return this->v_[X] + this->v_[Y] + this->v_[Z];
    }
    
    unsigned int index() const {
        // (y * (y + 1) + z * (2*y + z + 3)) /2 = ((y+z)*(y+z+1)+2z)/2
        const unsigned int yz = this->v_[Y] + this->v_[Z];
        return (yz * (yz+1) / 2 + this->v_[Z]);
    }

    bool isExist() const {
        return ((v_[0] >= 0) && (v_[1] >= 0) && (v_[2] >= 0));
    }

    /// 指定された軸[x(=0), y(=1), z(=3)] の要素を返す
    char get(const int i) const {
        assert((0 <= i) && (i < 3));
        return this->v_[i];
    }

    TlAngularMomentumVector operator+=(const TlAngularMomentumVector& rhs) {
        this->v_[X] += rhs.v_[X];
        this->v_[Y] += rhs.v_[Y];
        this->v_[Z] += rhs.v_[Z];

        return *this;
    }

    TlAngularMomentumVector operator-=(const TlAngularMomentumVector& rhs) {
        this->v_[X] -= rhs.v_[X];
        this->v_[Y] -= rhs.v_[Y];
        this->v_[Z] -= rhs.v_[Z];

        return *this;
    }
    
    friend TlAngularMomentumVector operator+(const TlAngularMomentumVector& rhs1,
                                             const TlAngularMomentumVector& rhs2) {
        TlAngularMomentumVector ans = rhs1;
        ans += rhs2;
        return ans;
    }

    friend TlAngularMomentumVector operator-(const TlAngularMomentumVector& rhs1,
                                             const TlAngularMomentumVector& rhs2) {
        TlAngularMomentumVector ans = rhs1;
        ans -= rhs2;
        return ans;
    }

    std::string debugOut() const {
        return TlUtils::format("[%d %d %d]", v_[0], v_[1], v_[2]);
    }
    
private:
    char v_[3];
};

#endif // TLANGULARMOMENTUMVECTOR_H
