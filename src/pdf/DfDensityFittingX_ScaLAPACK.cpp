#include "DfDensityFittingX_ScaLAPACK.h"

// SCALAPACK ===========================================================
DfDensityFittingX_ScaLAPACK::DfDensityFittingX_ScaLAPACK(TlSerializeData* pPdfParam)
    : DfDensityFittingTmpl<TlDistributeSymmetricMatrix, TlDistributeVector, DfEriX_Parallel>(pPdfParam)
{
}


DfDensityFittingX_ScaLAPACK::~DfDensityFittingX_ScaLAPACK()
{
}


void DfDensityFittingX_ScaLAPACK::exec()
{
    this->calc();
}

