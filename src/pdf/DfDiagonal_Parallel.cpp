#ifdef HAVE_CONFIG_H
#include "config.h"    // this file created by autotools
#endif // HAVE_CONFIG_H

#include "DfDiagonal_Parallel.h"
#include "CnError.h"
#include "TlCommunicate.h"
#include "TlVector.h"
#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"
#include "TlDistributeMatrix.h"
#include "TlDistributeSymmetricMatrix.h"
#include "TlUtils.h"

DfDiagonal_Parallel::DfDiagonal_Parallel(TlSerializeData* pPdfParam)
    : DfDiagonal(pPdfParam)
{
}

DfDiagonal_Parallel::~DfDiagonal_Parallel()
{
}

void DfDiagonal_Parallel::logger(const std::string& str) const
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfDiagonal::logger(str);
    }
}

void DfDiagonal_Parallel::DfDiagMain()
{
#ifdef HAVE_SCALAPACK
    if (this->m_bUsingSCALAPACK == true) {
        this->logger(" diagonal(parallel) using SCALAPACK.\n");
        this->DfDiagMain_SCALAPACK();
        return;
    }
#endif // HAVE_SCALAPACK

    // LAPACK 
    this->logger(" diagonal(parallel) using LAPACK.\n");
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfDiagonal::DfDiagMain();
    }
    rComm.barrier();
}

void DfDiagonal_Parallel::DfDiagQclo(const DfObject::RUN_TYPE runType, const std::string& fragname, int norbcut)
{
#ifdef HAVE_SCALAPACK
    if (this->m_bUsingSCALAPACK == true) {
        this->logger("DfDiagonal_Parallel using SCALAPACK.\n");
        this->DfDiagQclo_SCALAPACK(runType, fragname, norbcut);
        return;
    } else {
        this->logger("DfDiagonal_Parallel using LAPACK.\n");
    }
#endif // HAVE_SCALAPACK

    TlCommunicate& rComm = TlCommunicate::getInstance();

    if (rComm.isMaster() == true) {
        DfDiagonal::DfDiagQclo(runType, fragname, norbcut);
    }
    rComm.barrier();
}

void DfDiagonal_Parallel::DfDiagMain_SCALAPACK()
{
    switch(this->m_nMethodType) {
    case METHOD_RKS:
        this->main<TlDistributeMatrix, TlDistributeSymmetricMatrix>(RUN_RKS);
        break;

    case METHOD_UKS:
        this->main<TlDistributeMatrix, TlDistributeSymmetricMatrix>(RUN_UKS_ALPHA);
        this->main<TlDistributeMatrix, TlDistributeSymmetricMatrix>(RUN_UKS_BETA);
        break;

    case METHOD_ROKS:
        this->main<TlDistributeMatrix, TlDistributeSymmetricMatrix>(RUN_ROKS);
        break;

    default:
        CnErr.abort("DfDiagonal", "", "DfDiagMain", "the value of scftype is illegal");
        break;
    }
}

void DfDiagonal_Parallel::DfDiagQclo_SCALAPACK(const DfObject::RUN_TYPE runType,
                                               const std::string& fragname, int norbcut)
{
    this->m_nNumOfMOs = norbcut;
    this->main<TlDistributeMatrix, TlDistributeSymmetricMatrix>(runType, fragname, true);
}
