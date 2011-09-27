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
    int size() const {
        return (angularMomentum_+1)*(angularMomentum_+2) / 2;
    }

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

#endif // TLANGULARMOMENTUMVECTORSET_H
