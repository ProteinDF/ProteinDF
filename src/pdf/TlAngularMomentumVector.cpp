#include <iostream>
#include <cassert>

#include "TlAngularMomentumVector.h"
#include "TlUtils.h"

TlAngularMomentumVector::TlAngularMomentumVector(char x, char y, char z) {
    this->v_[X] = x;
    this->v_[Y] = y;
    this->v_[Z] = z;
}

TlAngularMomentumVector::TlAngularMomentumVector(const TlAngularMomentumVector& rhs) {
    this->v_[X] = rhs.v_[X];
    this->v_[Y] = rhs.v_[Y];
    this->v_[Z] = rhs.v_[Z];
}

int TlAngularMomentumVector::angularMomentum() const {
    return this->v_[X] + this->v_[Y] + this->v_[Z];
}

bool TlAngularMomentumVector::isExist() const {
    return ((v_[0] >= 0) && (v_[1] >= 0) && (v_[2] >= 0));
}

int TlAngularMomentumVector::get(const int i) const {
    assert((0 <= i) && (i < 3));
    return this->v_[i];
}

TlAngularMomentumVector TlAngularMomentumVector::operator+=(const TlAngularMomentumVector& rhs) {
    this->v_[X] += rhs.v_[X];
    this->v_[Y] += rhs.v_[Y];
    this->v_[Z] += rhs.v_[Z];
    
    return *this;
}

TlAngularMomentumVector TlAngularMomentumVector::operator-=(const TlAngularMomentumVector& rhs) {
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

