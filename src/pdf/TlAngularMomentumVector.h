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
                            char z = 0);

    TlAngularMomentumVector(const TlAngularMomentumVector& rhs);

    /// angular momentumを返す
    int angularMomentum() const;
    
    int index() const;

    bool isExist() const;

    /// 指定された軸[x(=0), y(=1), z(=3)] の要素を返す
    int get(const int i) const;

    TlAngularMomentumVector operator+=(const TlAngularMomentumVector& rhs);

    TlAngularMomentumVector operator-=(const TlAngularMomentumVector& rhs);
    
    friend TlAngularMomentumVector operator+(const TlAngularMomentumVector& rhs1,
                                             const TlAngularMomentumVector& rhs2);

    friend TlAngularMomentumVector operator-(const TlAngularMomentumVector& rhs1,
                                             const TlAngularMomentumVector& rhs2);

    std::string debugOut() const;

private:
    int v_[3];
};

#endif // TLANGULARMOMENTUMVECTOR_H
