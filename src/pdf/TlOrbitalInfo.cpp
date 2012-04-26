#include <iostream>
#include "TlOrbitalInfo.h"

TlOrbitalInfo::TlOrbitalInfo(const TlSerializeData& geomData,
                             const TlSerializeData& basisData)
{
    this->setCGTO(Fl_Gto(basisData));
    this->setAtoms(Fl_Geometry(geomData));
    this->makeOrbitalTable();
}


TlOrbitalInfo::~TlOrbitalInfo()
{
}

