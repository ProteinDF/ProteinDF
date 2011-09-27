#ifndef TLMOFIELD_H
#define TLMOFIELD_H

#include <vector>
#include "TlPosition.h"
#include "TlVector.h"
#include "TlSerializeData.h"
#include "TlOrbitalInfo.h"

class TlFieldData_Uniform;

class TlMoField {
public:
    TlMoField(const TlSerializeData& param);
    ~TlMoField();

    std::vector<double> makeMoFld(const TlVector& MO,
                                  const std::vector<TlPosition>& grids);

protected:
    double getPreFactor(const int nType, const TlPosition& pos);

protected:
    static const double INV_SQRT12;
    TlSerializeData param_;
    TlOrbitalInfo orbInfo_;
};

#endif // TLMOFIELD_H




