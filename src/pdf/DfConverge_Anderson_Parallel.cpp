#include "TlCommunicate.h"
#include "DfConverge_Anderson_Parallel.h"
#include "TlDistributeSymmetricMatrix.h"
#include "TlDistributeVector.h"

DfConverge_Anderson_Parallel::DfConverge_Anderson_Parallel(TlSerializeData* pPdfParam)
    : DfConverge_Anderson(pPdfParam)
{
}


DfConverge_Anderson_Parallel::~DfConverge_Anderson_Parallel()
{
}


void DfConverge_Anderson_Parallel::convergeRhoTilde()
{
#ifdef HAVE_SCALAPACK
    if (this->m_bUsingSCALAPACK == true) {
        this->log_.info("converge rho~ using ScaLAPACK.");
        this->convergeRhoTilde_ScaLAPACK();
        return;
    }
#endif // HAVE_SCALAPACK    
    this->log_.info("converge rho~ using LAPACK.");
    this->convergeRhoTilde_LAPACK();
}


void DfConverge_Anderson_Parallel::convergeRhoTilde_LAPACK()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfConverge_Anderson::convergeRhoTilde();
    }
}


void DfConverge_Anderson_Parallel::convergeRhoTilde_ScaLAPACK()
{
    switch (this->m_nMethodType) {
    case METHOD_RKS:
        DfConverge_Anderson::convergeRhoTilde<TlDistributeVector>(DfObject::RUN_RKS);
        break;
    case METHOD_UKS:
        DfConverge_Anderson::convergeRhoTilde<TlDistributeVector>(DfObject::RUN_UKS_ALPHA);
        DfConverge_Anderson::convergeRhoTilde<TlDistributeVector>(DfObject::RUN_UKS_BETA);
        break;
    case METHOD_ROKS:
        DfConverge_Anderson::convergeRhoTilde<TlDistributeVector>(DfObject::RUN_UKS_ALPHA);
        DfConverge_Anderson::convergeRhoTilde<TlDistributeVector>(DfObject::RUN_UKS_BETA);
        break;
    default:
        std::cerr << "program error. @DfConverge_Anderson_Parallel::convergeRhoTilde()" << std::endl;
        break;
    }
}


void DfConverge_Anderson_Parallel::convergeKSMatrix()
{
#ifdef HAVE_SCALAPACK
    if (this->m_bUsingSCALAPACK == true) {
        this->log_.info("converge KS matrix using ScaLAPACK.");
        this->convergeKSMatrix_ScaLAPACK();
        return;
    }
#endif // HAVE_SCALAPACK    
    this->log_.info("converge KS matrix using LAPACK.");
    this->convergeKSMatrix_LAPACK();
}


void DfConverge_Anderson_Parallel::convergeKSMatrix_LAPACK()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfConverge_Anderson::convergeKSMatrix();
    }
}


void DfConverge_Anderson_Parallel::convergeKSMatrix_ScaLAPACK()
{
    switch (this->m_nMethodType) {
    case METHOD_RKS:
        DfConverge_Anderson::convergeKSMatrix<TlDistributeSymmetricMatrix, TlDistributeVector>(DfObject::RUN_RKS);
        break;
    case METHOD_UKS:
        DfConverge_Anderson::convergeKSMatrix<TlDistributeSymmetricMatrix, TlDistributeVector>(DfObject::RUN_UKS_ALPHA);
        DfConverge_Anderson::convergeKSMatrix<TlDistributeSymmetricMatrix, TlDistributeVector>(DfObject::RUN_UKS_BETA);
        break;
    case METHOD_ROKS:
        DfConverge_Anderson::convergeKSMatrix<TlDistributeSymmetricMatrix, TlDistributeVector>(DfObject::RUN_UKS_ALPHA);
        DfConverge_Anderson::convergeKSMatrix<TlDistributeSymmetricMatrix, TlDistributeVector>(DfObject::RUN_UKS_BETA);
        break;
    default:
        std::cerr << "program error. @DfConverge_Anderson_Parallel::convergeKSMatrix()" << std::endl;
        break;
    }
}


void DfConverge_Anderson_Parallel::convergePMatrix()
{
#ifdef HAVE_SCALAPACK
    if (this->m_bUsingSCALAPACK == true) {
        this->log_.info("converge density matrix using ScaLAPACK.");
        this->convergePMatrix_ScaLAPACK();
        return;
    }
#endif // HAVE_SCALAPACK    
    this->log_.info("converge density matrix using LAPACK.");
    this->convergePMatrix_LAPACK();
}


void DfConverge_Anderson_Parallel::convergePMatrix_LAPACK()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfConverge_Anderson::convergePMatrix();
    }
}


void DfConverge_Anderson_Parallel::convergePMatrix_ScaLAPACK()
{
    switch (this->m_nMethodType) {
    case METHOD_RKS:
        DfConverge_Anderson::convergePMatrix<TlDistributeSymmetricMatrix, TlDistributeVector>(DfObject::RUN_RKS);
        break;
    case METHOD_UKS:
        DfConverge_Anderson::convergePMatrix<TlDistributeSymmetricMatrix, TlDistributeVector>(DfObject::RUN_UKS_ALPHA);
        DfConverge_Anderson::convergePMatrix<TlDistributeSymmetricMatrix, TlDistributeVector>(DfObject::RUN_UKS_BETA);
        break;
    case METHOD_ROKS:
        DfConverge_Anderson::convergePMatrix<TlDistributeSymmetricMatrix, TlDistributeVector>(DfObject::RUN_UKS_ALPHA);
        DfConverge_Anderson::convergePMatrix<TlDistributeSymmetricMatrix, TlDistributeVector>(DfObject::RUN_UKS_BETA);
        break;
    default:
        std::cerr << "program error. @DfConverge_Anderson_Distribute::convergePMatrix()" << std::endl;
        break;
    }
}

