#include <iostream>
#include "TlOrbitalInfo.h"

#include "Fl_Geometry.h"
#include "Fl_Gto_Orbital.h"

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

