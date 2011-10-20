#ifdef HAVE_CONFIG_H
#include "config.h"    // this file created by autotools
#endif // HAVE_CONFIG_H

#include "DfDensityFittingX.h"

DfDensityFittingX::DfDensityFittingX(TlSerializeData* pPdfParam)
    : DfDensityFittingTmpl<TlSymmetricMatrix, TlVector, DfEriX>(pPdfParam)
{
}


DfDensityFittingX::~DfDensityFittingX()
{
}


void DfDensityFittingX::exec()
{
    this->calc();
}

