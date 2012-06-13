#include "DfInitialGuessHarris_Parallel.h"
#include "TlCommunicate.h"
#include "DfOverlapX_Parallel.h"
#include "DfPopulation_Parallel.h"

DfInitialGuessHarris_Parallel::DfInitialGuessHarris_Parallel(TlSerializeData* pPdfParam)
    : DfInitialGuessHarris(pPdfParam)
{
}


DfInitialGuessHarris_Parallel::~DfInitialGuessHarris_Parallel()
{
}


void DfInitialGuessHarris_Parallel::main()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
#ifdef HAVE_SCALAPACK
    if (this->m_bUsingSCALAPACK) {
        DfInitialGuessHarris::calcInitialDensityMatrix<TlDistributeMatrix,
            TlDistributeSymmetricMatrix,
            DfOverlapX_Parallel, 
            DfPopulation_Parallel>();
    } else {
        if (rComm.isMaster()) {
            DfInitialGuessHarris::main();
        }
    }
#else
    if (rComm.isMaster()) {
        DfInitialGuessHarris::main();
    }
#endif // HAVE_SCALAPACK
}
