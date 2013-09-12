#include "DfGridFreeXC_Parallel.h"
#include "DfOverlapX_Parallel.h"
#include "DfXMatrix_Parallel.h"

DfGridFreeXC_Parallel::DfGridFreeXC_Parallel(TlSerializeData* pPdfParam)
    : DfGridFreeXC(pPdfParam) {
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

