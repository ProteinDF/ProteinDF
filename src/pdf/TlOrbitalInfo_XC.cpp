#include <iostream>
#include "TlOrbitalInfo_XC.h"

#include "Fl_Geometry.h"

TlOrbitalInfo_XC::TlOrbitalInfo_XC(const TlSerializeData& geomData,
                                   const TlSerializeData& basisData)
    : TlOrbitalInfoObject() {
    this->setCGTO(Fl_Gto(basisData));
    this->setAtoms(Fl_Geometry(geomData));
    this->makeOrbitalTable();
}


TlOrbitalInfo_XC::~TlOrbitalInfo_XC()
{
}

