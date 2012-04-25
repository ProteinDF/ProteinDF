#include "TlCommunicate.h"
#include "DfDiffDensityMatrix_Parallel.h"
#include "TlDistributeSymmetricMatrix.h"
#include "TlFile.h"

DfDiffDensityMatrix_Parallel::DfDiffDensityMatrix_Parallel(TlSerializeData* pPdfParam)
    : DfDiffDensityMatrix(pPdfParam)
{
}


DfDiffDensityMatrix_Parallel::~DfDiffDensityMatrix_Parallel()
{
}


void DfDiffDensityMatrix_Parallel::exec()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

#ifdef HAVE_SCALAPACK
    if (this->m_bUsingSCALAPACK == true) {
        // using ScaLAPACK
        switch (this->m_nMethodType) {
        case METHOD_RKS:
            this->calc_usingScaLAPACK(RUN_RKS, this->m_nIteration);
            break;

        case METHOD_UKS:
            this->calc_usingScaLAPACK(RUN_UKS_ALPHA, this->m_nIteration);
            this->calc_usingScaLAPACK(RUN_UKS_BETA, this->m_nIteration);
            break;

        case METHOD_ROKS:
            DfDiffDensityMatrix::calc_ROKS<TlDistributeSymmetricMatrix>();
            break;

        default:
            std::cerr << "program error. " << __FILE__ << ":" << __LINE__ << std::endl;
            abort();
        }

    } else {
        // using LAPACK
        if (rComm.isMaster() == true) {
            DfDiffDensityMatrix::exec();
        }
        rComm.barrier();
    }
#else
    {
        // using LAPACK
        if (rComm.isMaster() == true) {
            DfDiffDensityMatrix::exec();
        }
        rComm.barrier();
    }
#endif // HAVE_SCALAPACK
}


void DfDiffDensityMatrix_Parallel::calc_usingScaLAPACK(const DfObject::RUN_TYPE runType,
                                                       const int iteration)
{
    this->log_.info("(delta P) is build using ScaLAPACK.");
    TlDistributeSymmetricMatrix P = this->getPpqMatrix<TlDistributeSymmetricMatrix>(runType, iteration -1);
    if (TlFile::isExist(this->getPpqMatrixPath(runType, iteration -2)) == true) {
        P -= (this->getPpqMatrix<TlDistributeSymmetricMatrix>(runType, iteration -2));
    }

    this->saveDiffDensityMatrix(runType, iteration, P);
}
