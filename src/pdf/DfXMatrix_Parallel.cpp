#ifdef HAVE_CONFIG_H
#include "config.h"    // this file created by autotools
#endif // HAVE_CONFIG_H

#include "DfXMatrix_Parallel.h"
#include "TlCommunicate.h"
#include "TlDistributeMatrix.h"
#include "TlDistributeSymmetricMatrix.h"

#define FAST_TRANCATE

DfXMatrix_Parallel::DfXMatrix_Parallel(TlSerializeData* pPdfParam)
    : DfXMatrix(pPdfParam)
{
}

DfXMatrix_Parallel::~DfXMatrix_Parallel()
{
}

void DfXMatrix_Parallel::buildX()
{
#ifdef HAVE_SCALAPACK
    if (this->m_bUsingSCALAPACK == true) {
        this->buildX_ScaLAPACK();
    } else {
        this->buildX_LAPACK();
    }
#else
    {
        this->buildX_LAPACK();
    }
#endif // HAVE_SCALAPACK
}

void DfXMatrix_Parallel::buildX_LAPACK()
{
    this->log_.info("build X using replica mode.");

    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfXMatrix::buildX();
    }

    index_type numOfMOs = 0;
    if (rComm.isMaster() == true) {
        numOfMOs = (*(this->pPdfParam_))["num_of_MOs"].getInt();
    }
    rComm.broadcast(numOfMOs);
    (*(this->pPdfParam_))["num_of_MOs"] = numOfMOs;
}

void DfXMatrix_Parallel::buildX_ScaLAPACK()
{
    TlDistributeSymmetricMatrix S = this->getSpqMatrix<TlDistributeSymmetricMatrix>();
    TlDistributeMatrix X;
    TlDistributeMatrix Xinv;
    
    this->canonicalOrthogonalize(S, &X, &Xinv);

    (*(this->pPdfParam_))["num_of_MOs"] = X.getNumOfCols();

    DfObject::saveXMatrix(X);
    DfObject::saveXInvMatrix(Xinv);
}




