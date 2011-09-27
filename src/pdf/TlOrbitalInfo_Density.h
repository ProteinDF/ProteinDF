#ifndef TLORBITAL_INFO_DENSITY_H
#define TLORBITAL_INFO_DENSITY_H

#include "TlOrbitalInfoObject.h"
#include "Fl_Gto_Density.h"

class TlOrbitalInfo_Density : public TlOrbitalInfoObject {
public:
    TlOrbitalInfo_Density(const TlSerializeData& geomData,
                          const TlSerializeData& basisData);

    virtual ~TlOrbitalInfo_Density();
};

#endif // TLORBITAL_INFO_DENSITY_H
