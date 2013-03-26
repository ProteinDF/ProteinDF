#ifdef HAVE_CONFIG_H
#include "config.h"    // this file created by autotools
#endif // HAVE_CONFIG_H

#include "DfInvMatrix_Parallel.h"
#include "TlCommunicate.h"
#include "TlSymmetricMatrix.h"
#include "TlDistributeSymmetricMatrix.h"

DfInvMatrix_Parallel::DfInvMatrix_Parallel(TlSerializeData* pPdfParam)
    : DfInvMatrix(pPdfParam)
{
}


DfInvMatrix_Parallel::~DfInvMatrix_Parallel()
{
}


void DfInvMatrix_Parallel::DfInvMain()
{
#ifdef HAVE_SCALAPACK
    if (this->m_bUsingSCALAPACK == true) {
        this->exec<TlDistributeSymmetricMatrix>();
        return;
    }
#endif // HAVE_SCALAPACK
    
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfInvMatrix::DfInvMain();
    }
    rComm.barrier();
}




