#include "DfConverge_Damping_Parallel.h"
#include "TlDistributeVector.h"
#include "TlDistributeSymmetricMatrix.h"

DfConverge_Damping_Parallel::DfConverge_Damping_Parallel(TlSerializeData* pPdfParam)
    : DfConverge_Damping(pPdfParam)
{
}


DfConverge_Damping_Parallel::~DfConverge_Damping_Parallel()
{
}


void DfConverge_Damping_Parallel::convergeRhoTilde()
{
#ifdef HAVE_SCALAPACK
    if (this->m_bUsingSCALAPACK == true) {
        this->log_.info(" converge rho~ using ScaLAPACK.");
        this->convergeRhoTilde_ScaLAPACK();
        return;
    }
#endif // HAVE_SCALAPACK    
    this->log_.info(" converge rho~ using LAPACK.");
    this->convergeRhoTilde_LAPACK();
}


void DfConverge_Damping_Parallel::convergeRhoTilde_LAPACK()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfConverge_Damping::convergeRhoTilde();
    }
}


void DfConverge_Damping_Parallel::convergeRhoTilde_ScaLAPACK()
{
    switch (this->m_nMethodType) {
    case METHOD_RKS:
        DfConverge_Damping::convergeRhoTilde<TlDistributeVector>(DfObject::RUN_RKS);
        break;
    case METHOD_UKS:
        DfConverge_Damping::convergeRhoTilde<TlDistributeVector>(DfObject::RUN_UKS_ALPHA);
        DfConverge_Damping::convergeRhoTilde<TlDistributeVector>(DfObject::RUN_UKS_BETA);
        break;
    case METHOD_ROKS:
        DfConverge_Damping::convergeRhoTilde<TlDistributeVector>(DfObject::RUN_UKS_ALPHA);
        DfConverge_Damping::convergeRhoTilde<TlDistributeVector>(DfObject::RUN_UKS_BETA);
        break;
    default:
        std::cerr << "program error. @DfConverge_Damping_Parallel::convergeRhoTilde()" << std::endl;
        break;
    }
}


void DfConverge_Damping_Parallel::convergeKSMatrix()
{
#ifdef HAVE_SCALAPACK
    if (this->m_bUsingSCALAPACK == true) {
        this->log_.info(" converge KS matrix using ScaLAPACK.");
        this->convergeKSMatrix_ScaLAPACK();
        return;
    }
#endif // HAVE_SCALAPACK
    this->log_.info(" converge KS matrix using LAPACK.");
    this->convergeKSMatrix_LAPACK();
}


void DfConverge_Damping_Parallel::convergeKSMatrix_LAPACK()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfConverge_Damping::convergeKSMatrix();
    }
}


void DfConverge_Damping_Parallel::convergeKSMatrix_ScaLAPACK()
{
    switch (this->m_nMethodType) {
    case METHOD_RKS:
        DfConverge_Damping::convergeKSMatrix<TlDistributeSymmetricMatrix>(DfObject::RUN_RKS);
        break;
    case METHOD_UKS:
        DfConverge_Damping::convergeKSMatrix<TlDistributeSymmetricMatrix>(DfObject::RUN_UKS_ALPHA);
        DfConverge_Damping::convergeKSMatrix<TlDistributeSymmetricMatrix>(DfObject::RUN_UKS_BETA);
        break;
    case METHOD_ROKS:
        DfConverge_Damping::convergeKSMatrix<TlDistributeSymmetricMatrix>(DfObject::RUN_UKS_ALPHA);
        DfConverge_Damping::convergeKSMatrix<TlDistributeSymmetricMatrix>(DfObject::RUN_UKS_BETA);
        break;
    default:
        std::cerr << "program error. @DfConverge_Damping_Parallel::convergeKSMatrix()" << std::endl;
        break;
    }
}


void DfConverge_Damping_Parallel::convergePMatrix()
{
#ifdef HAVE_SCALAPACK
    if (this->m_bUsingSCALAPACK == true) {
        this->log_.info(" converge density matrix using ScaLAPACK.");
        this->convergePMatrix_ScaLAPACK();
        return;
    }
#endif // HAVE_SCALAPACK
    this->log_.info(" converge density matrix using LAPACK.");
    this->convergePMatrix_LAPACK();
}


void DfConverge_Damping_Parallel::convergePMatrix_LAPACK()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfConverge_Damping::convergePMatrix();
    }
    rComm.barrier();
}


void DfConverge_Damping_Parallel::convergePMatrix_ScaLAPACK()
{
    switch (this->m_nMethodType) {
    case METHOD_RKS:
        DfConverge_Damping::convergePMatrix<TlDistributeSymmetricMatrix>(DfObject::RUN_RKS);
        break;
    case METHOD_UKS:
        DfConverge_Damping::convergePMatrix<TlDistributeSymmetricMatrix>(DfObject::RUN_UKS_ALPHA);
        DfConverge_Damping::convergePMatrix<TlDistributeSymmetricMatrix>(DfObject::RUN_UKS_BETA);
        break;
    case METHOD_ROKS:
        DfConverge_Damping::convergePMatrix<TlDistributeSymmetricMatrix>(DfObject::RUN_UKS_ALPHA);
        DfConverge_Damping::convergePMatrix<TlDistributeSymmetricMatrix>(DfObject::RUN_UKS_BETA);
        break;
    default:
        std::cerr << "program error. @DfConverge_Damping_Distribute::convergePMatrix()" << std::endl;
        break;
    }
}

