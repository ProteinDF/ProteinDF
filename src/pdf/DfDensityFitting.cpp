#ifdef HAVE_CONFIG_H
#include "config.h"    // this file created by autotools
#endif // HAVE_CONFIG_H

#include "DfDensityFitting.h"

DfDensityFitting::DfDensityFitting(TlSerializeData* pPdfParam)
    : DfDensityFittingTmpl<TlSymmetricMatrix, TlVector, DfEri>(pPdfParam)
{
}


DfDensityFitting::~DfDensityFitting()
{
}


void DfDensityFitting::exec()
{
    this->calc();
}


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

