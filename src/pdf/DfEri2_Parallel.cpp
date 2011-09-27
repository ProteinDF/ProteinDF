#include "DfEri2_Parallel.h"
#include "DfEri_Parallel.h"
#include "TlCommunicate.h"

DfEri2_Parallel::DfEri2_Parallel(TlSerializeData* pPdfParam)
        : DfEri2(pPdfParam)
{
}

DfEri2_Parallel::~DfEri2_Parallel()
{
}

DfEri* DfEri2_Parallel::getDfEri()
{
    DfEri* pDfEri = new DfEri_Parallel(this->pPdfParam_);

    return pDfEri;
}

TlSymmetricMatrix DfEri2_Parallel::loadInvSquareVMatrix()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    TlSymmetricMatrix V;
    if (rComm.isMaster() == true) {
        V = DfEri2::loadInvSquareVMatrix();
    }

    rComm.broadcast(V);

    return V;
}

TlMatrix DfEri2_Parallel::loadCutoffCMatrix()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    TlMatrix C;

    if (rComm.isMaster() == true) {
        C = DfEri2::loadCutoffCMatrix();
    }

    rComm.broadcast(C);

    return C;
}

TlSymmetricMatrix DfEri2_Parallel::getPMatrix(const int iteration)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    TlSymmetricMatrix P;
    if (rComm.isMaster() == true) {
        P = DfEri2::getPMatrix(iteration);
    }

    rComm.broadcast(P);

    return P;
}

TlSymmetricMatrix DfEri2_Parallel::loadKMatrix(const int nIteration)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    TlSymmetricMatrix K;
    if (rComm.isMaster() == true) {
        K = DfEri2::loadKMatrix(nIteration);
    }

    rComm.broadcast(K);

    return K;
}

void DfEri2_Parallel::saveKMatrix(const TlSymmetricMatrix& K, const int nIteration)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    if (rComm.isMaster() == true) {
        DfEri2::saveKMatrix(K, nIteration);
    }

    rComm.barrier();
}

