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


void DfDensityFittingX_ScaLAPACK::logger(const std::string& str) const
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    if (rComm.isMaster() == true) {
        DfObject::logger(str);
    }
}


