#ifndef TLORBITALINFO_H
#define TLORBITALINFO_H

#include "TlOrbitalInfoObject.h"

class TlOrbitalInfo : public TlOrbitalInfoObject {
public:
    TlOrbitalInfo(const TlSerializeData& geomData,
                  const TlSerializeData& basisData);

    virtual ~TlOrbitalInfo();
};

#endif // TLORBITALINFO_H
