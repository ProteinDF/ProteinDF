#include <cassert>
#include "TlAngularMomentumVectorSet.h"

const TlAngularMomentumVectorSet::AMV_DB TlAngularMomentumVectorSet::db_;

TlAngularMomentumVectorSet::TlAngularMomentumVectorSet(const int angularMomentum)
    : angularMomentum_(angularMomentum) {
    assert((0 <= angularMomentum) && (angularMomentum < MAX_ANGULAR_MOMENTUM));
}


// int TlAngularMomentumVectorSet::size() const
// {
//     return this->db_.data[this->angularMomentum_].size();
// }


const TlAngularMomentumVector& TlAngularMomentumVectorSet::get(const int index) const
{
    assert((0 <= index) && (index < this->size()));
    return this->db_.data[this->angularMomentum_][index];
}

