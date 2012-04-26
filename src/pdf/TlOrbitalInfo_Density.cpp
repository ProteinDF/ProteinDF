#include <iostream>
#include "TlOrbitalInfo_Density.h"

TlOrbitalInfo_Density::TlOrbitalInfo_Density(const TlSerializeData& geomData,
                                             const TlSerializeData& basisData)
{
    this->setCGTO_coulomb(Fl_Gto(basisData));
    this->setAtoms(Fl_Geometry(geomData));
    this->makeOrbitalTable();
}


TlOrbitalInfo_Density::~TlOrbitalInfo_Density()
{
}

