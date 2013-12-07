#include "DfGridFreeXC_Parallel.h"
#include "DfOverlapX_Parallel.h"
#include "DfXMatrix_Parallel.h"
#include "DfCD_Parallel.h"

DfGridFreeXC_Parallel::DfGridFreeXC_Parallel(TlSerializeData* pPdfParam)
    : DfGridFreeXC(pPdfParam) {
    this->log_.info("DfGridFreeXC_Parallel::DfGridFreeXC_Parallel()");
}

DfGridFreeXC_Parallel::~DfGridFreeXC_Parallel()
{
}

DfOverlapX* DfGridFreeXC_Parallel::getDfOverlapObject()
{
    DfOverlapX* pDfOverlapX = new DfOverlapX_Parallel(this->pPdfParam_);
    return pDfOverlapX;
}

DfXMatrix* DfGridFreeXC_Parallel::getDfXMatrixObject()
{
    DfXMatrix* pDfXMatrix = new DfXMatrix_Parallel(this->pPdfParam_);
    return pDfXMatrix;
}

// before SCF ==================================================================
void DfGridFreeXC_Parallel::preprocessBeforeSCF()
{
#ifdef HAVE_SCALAPACK
    if (this->m_bUsingSCALAPACK == true) {
        this->preprocessBeforeSCF_ScaLAPACK();
    } else {
        this->preprocessBeforeSCF_LAPACK();
    }
#else
    {
        this->preprocessBeforeSCF_LAPACK();
    }
#endif // HAVE_SCALAPACK
}

void DfGridFreeXC_Parallel::preprocessBeforeSCF_LAPACK()
{
    this->preprocessBeforeSCF_templ<DfOverlapX_Parallel,
                                    DfXMatrix_Parallel,
                                    TlSymmetricMatrix, TlMatrix>();
}

void DfGridFreeXC_Parallel::preprocessBeforeSCF_ScaLAPACK()
{
    this->preprocessBeforeSCF_templ<DfOverlapX_Parallel,
                                    DfXMatrix_Parallel,
                                    TlDistributeSymmetricMatrix,
                                    TlDistributeMatrix>();
}


// in SCF ======================================================================
void DfGridFreeXC_Parallel::buildFxc_LDA()
{
    this->log_.info("DfGridFreeXC_Parallel::buildFxc_LDA()");
#ifdef HAVE_SCALAPACK
    if (this->m_bUsingSCALAPACK == true) {
        this->buildFxc_LDA_ScaLAPACK();
    } else {
        this->buildFxc_LDA_LAPACK();
    }
#else
    {
        this->buildFxc_LDA_LAPACK();
    }
#endif // HAVE_SCALAPACK
}

void DfGridFreeXC_Parallel::buildFxc_LDA_LAPACK()
{
    this->log_.info("DfGridFreeXC_Parallel::buildFxc_LDA_LAPACK()");
    DfGridFreeXC::buildFxc_LDA_method<DfOverlapX_Parallel,
                                      DfCD_Parallel,
                                      TlSymmetricMatrix,
                                      TlMatrix>();
}

void DfGridFreeXC_Parallel::buildFxc_LDA_ScaLAPACK()
{
    this->log_.info("DfGridFreeXC_Parallel::buildFxc_LDA_ScaLAPACK()");
    DfGridFreeXC::buildFxc_LDA_method<DfOverlapX_Parallel,
                                      DfCD_Parallel,
                                      TlDistributeSymmetricMatrix,
                                      TlDistributeMatrix>();
}


