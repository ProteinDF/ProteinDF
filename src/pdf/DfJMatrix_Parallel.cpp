#include <cassert>
#include "TlCommunicate.h"
#include "DfJMatrix_Parallel.h"
#include "DfEriX_Parallel.h"

DfJMatrix_Parallel::DfJMatrix_Parallel(TlSerializeData* pPdfParam)
    : DfJMatrix(pPdfParam) {
}


DfJMatrix_Parallel::~DfJMatrix_Parallel()
{
}


void DfJMatrix_Parallel::buildJMatrix()
{
    const int iteration = this->m_nIteration;
    const RUN_TYPE runType = RUN_RKS;

#ifdef HAVE_SCALAPACK
    if (this->m_bUsingSCALAPACK == true) {
        this->buildJMatrix_ScaLAPACK(iteration, runType);
    } else {
        this->buildJMatrix_LAPACK(iteration, runType);
    }
#else
    this->buildJMatrix_LAPACK(iteration, runType);
#endif // HAVE_SCALAPACK
}


void DfJMatrix_Parallel::buildJMatrix_LAPACK(const int iteration,
                                             const RUN_TYPE runType)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    
    const index_type numOfAOs = this->m_nNumOfAOs;
    
    TlSymmetricMatrix P;
    if (rComm.isMaster() == true) {
        P = this->getDiffDensityMatrix<TlSymmetricMatrix>(runType, iteration);
    }
    rComm.broadcast(P);
    assert(P.getNumOfRows() == numOfAOs);
    
    DfEriX_Parallel dfEri(this->pPdfParam_);
    TlSymmetricMatrix J(numOfAOs);
    dfEri.getJpq(P, &J);

    if (rComm.isMaster() == true) {
        if (iteration >= 2) {
            TlSymmetricMatrix prevJ = this->getJMatrix<TlSymmetricMatrix>(iteration -1);
            J += prevJ;
        }

        this->saveJMatrix(iteration, J);
    }
    rComm.barrier();
}


void DfJMatrix_Parallel::buildJMatrix_ScaLAPACK(const int iteration,
                                                const RUN_TYPE runType)
{
    const index_type numOfAOs = this->m_nNumOfAOs;
    const TlDistributeSymmetricMatrix P =
        this->getDiffDensityMatrix<TlDistributeSymmetricMatrix>(runType, iteration);

    DfEriX_Parallel dfEri(this->pPdfParam_);
    TlDistributeSymmetricMatrix J(numOfAOs);
    dfEri.getJpq_D(P, &J);

    if (iteration >= 2) {
        TlDistributeSymmetricMatrix prevJ = this->getJMatrix<TlDistributeSymmetricMatrix>(iteration -1);
        J += prevJ;
    }

    this->saveJMatrix(iteration, J);
}

