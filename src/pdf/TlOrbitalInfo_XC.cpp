#include <iostream>
#include "TlOrbitalInfo_XC.h"

#include "Fl_Geometry.h"
#include "Fl_Gto_Xcpot.h"

TlOrbitalInfo_XC::TlOrbitalInfo_XC(const TlSerializeData& geomData,
                                   const TlSerializeData& basisData)
{
    this->setCGTO(Fl_Gto(basisData));
    this->setAtoms(Fl_Geometry(geomData));
    this->makeOrbitalTable();
}


TlOrbitalInfo_XC::~TlOrbitalInfo_XC()
{
}

