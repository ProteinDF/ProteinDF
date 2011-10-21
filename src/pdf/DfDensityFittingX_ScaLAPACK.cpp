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


template<>
TlDistributeVector
DfDensityFittingTmpl<TlDistributeSymmetricMatrix, TlDistributeVector, DfEriX_Parallel>::getTalpha(const RUN_TYPE runType)
{
    TlDistributeVector t_alpha;
    std::string suffix = "";
    if (runType == RUN_UKS_ALPHA) {
        suffix = "a";
    } else if (runType == RUN_UKS_BETA) {
        suffix = "b";
    }

    if (this->m_bDiskUtilization == false) {
        TlDistributeSymmetricMatrix P =
            DfObject::getDiffDensityMatrix<TlDistributeSymmetricMatrix>(runType, this->m_nIteration);
        t_alpha = this->calcTAlpha_DIRECT(P);
        t_alpha += this->getTalpha(runType, this->m_nIteration -1);
    } else {
        abort();
    }

    // save
    if (this->m_bDiskUtilization == false) {
        t_alpha.save("fl_Work/fl_Vct_Talpha" + suffix + TlUtils::xtos(this->m_nIteration));
    }

    return t_alpha;
}

