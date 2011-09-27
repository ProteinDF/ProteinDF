#include "DfPopulation_Parallel.h"
#include "TlDistributeSymmetricMatrix.h"
#include "TlCommunicate.h"

DfPopulation_Parallel::DfPopulation_Parallel(TlSerializeData* pPdfParam)
    : DfPopulation(pPdfParam)
{
}

DfPopulation_Parallel::~DfPopulation_Parallel()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    rComm.barrier();
}


void DfPopulation_Parallel::calcPop(const int iteration)
{
#ifdef HAVE_SCALAPACK
    if (this->m_bUsingSCALAPACK == true) {
        DfPopulation::calcPop<TlDistributeSymmetricMatrix>(iteration);
        return;
    }
#endif // HAVE_SCALAPACK

    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfPopulation::calcPop(iteration);
    }
    rComm.barrier();
}


