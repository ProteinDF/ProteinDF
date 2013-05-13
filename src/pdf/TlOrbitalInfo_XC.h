#ifndef TLORBITAL_INFO_XC_H
#define TLORBITAL_INFO_XC_H

#include "TlOrbitalInfoObject.h"

class TlOrbitalInfo_XC : public TlOrbitalInfoObject {
public:
    TlOrbitalInfo_XC(const TlSerializeData& geomData,
                     const TlSerializeData& basisData);

    virtual ~TlOrbitalInfo_XC();
};

#endif // TLORBITAL_INFO_XC_H
