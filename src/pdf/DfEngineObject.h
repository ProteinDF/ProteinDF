#ifndef DFENGINEOBJECT_H
#define DFENGINEOBJECT_H

#include <vector>
#include "TlOrbitalInfoObject.h"

// s=0, p=1, d=2
#define AM_LEVEL (2)
// 2次微分まで対応
#define AM_GRAD_LEVEL (2)

class DfEngineObject {
public:
    typedef int index_type;
    
public:
    void setPrimitiveLevelThreshold(const double threshold) {
        this->primitiveLevelThreshold_ = threshold;
    }
    double getPrimitiveLevelThreshold() const {
        return this->primitiveLevelThreshold_;
    }
    
    virtual void calc(const int diff1, const TlOrbitalInfoObject& orbInfo1, const index_type shell1,
                      const int diff2, const TlOrbitalInfoObject& orbInfo2, const index_type shell2,
                      const int diff3, const TlOrbitalInfoObject& orbInfo3, const index_type shell3,
                      const int diff4, const TlOrbitalInfoObject& orbInfo4, const index_type shell4) =0;

public:
    virtual ~DfEngineObject() {
    };

public:
    virtual double value(const index_type index) const =0;

protected:
    // see J. Chem. Phys., 105, 2726 (1996), eq33
    double primitiveLevelThreshold_;
};

#endif // DFENGINEOBJECT_H

