#ifndef TLORBITALINFO_H
#define TLORBITALINFO_H

#include "TlOrbitalInfoObject.h"

class TlOrbitalInfo : public TlOrbitalInfoObject {
public:
    //explicit TlOrbitalInfo(const std::string& baseDir = ".");
    TlOrbitalInfo(const TlSerializeData& geomData,
                  const TlSerializeData& basisData);
    // TlOrbitalInfo(const Fl_Geometry& flGeom,
    //               const Fl_Gto& flGto);

    virtual ~TlOrbitalInfo();
};

#endif // TLORBITALINFO_H
