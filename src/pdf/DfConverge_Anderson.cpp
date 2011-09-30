#include "DfConverge_Anderson.h"
#include "TlSymmetricMatrix.h"
#include "TlLogX.h"

DfConverge_Anderson::DfConverge_Anderson(TlSerializeData* pPdfParam)
    : DfConverge_Damping(pPdfParam)
{
    const TlSerializeData& pdfParam = *pPdfParam;
    this->m_nStartIterationOfAnderson =
        std::max(pdfParam["scf-acceleration/anderson/start-number"].getInt(), 3);

    this->m_dDampingFactorOfAnderson =
        pdfParam["scf-acceleration/anderson/damping-factor"].getDouble();

}

DfConverge_Anderson::~DfConverge_Anderson()
{
}


void DfConverge_Anderson::convergeRhoTilde()
{
    switch (this->m_nMethodType) {
    case METHOD_RKS:
        this->convergeRhoTilde<TlVector>(DfObject::RUN_RKS);
        break;
    case METHOD_UKS:
        this->convergeRhoTilde<TlVector>(DfObject::RUN_UKS_ALPHA);
        this->convergeRhoTilde<TlVector>(DfObject::RUN_UKS_BETA);
        break;
    case METHOD_ROKS:
        this->convergeRhoTilde<TlVector>(DfObject::RUN_UKS_ALPHA);
        this->convergeRhoTilde<TlVector>(DfObject::RUN_UKS_BETA);
        break;
    default:
        std::cerr << "program error. @DfConverge_Anderson::convergeRhoTilde()" << std::endl;
        break;
    }
}

void DfConverge_Anderson::convergeKSMatrix()
{
    switch (this->m_nMethodType) {
    case METHOD_RKS:
        this->convergeKSMatrix<TlSymmetricMatrix, TlVector>(DfObject::RUN_RKS);
        break;
    case METHOD_UKS:
        this->convergeKSMatrix<TlSymmetricMatrix, TlVector>(DfObject::RUN_UKS_ALPHA);
        this->convergeKSMatrix<TlSymmetricMatrix, TlVector>(DfObject::RUN_UKS_BETA);
        break;
    case METHOD_ROKS:
        this->convergeKSMatrix<TlSymmetricMatrix, TlVector>(DfObject::RUN_UKS_ALPHA);
        this->convergeKSMatrix<TlSymmetricMatrix, TlVector>(DfObject::RUN_UKS_BETA);
        break;
    default:
        std::cerr << "program error. @DfConverge_Anderson::convergeKSMatrix()" << std::endl;
        break;
    }
}

void DfConverge_Anderson::convergePMatrix()
{
    switch (this->m_nMethodType) {
    case METHOD_RKS:
        this->convergePMatrix<TlSymmetricMatrix, TlVector>(DfObject::RUN_RKS);
        break;
    case METHOD_UKS:
        this->convergePMatrix<TlSymmetricMatrix, TlVector>(DfObject::RUN_UKS_ALPHA);
        this->convergePMatrix<TlSymmetricMatrix, TlVector>(DfObject::RUN_UKS_BETA);
        break;
    case METHOD_ROKS:
        this->convergePMatrix<TlSymmetricMatrix, TlVector>(DfObject::RUN_UKS_ALPHA);
        this->convergePMatrix<TlSymmetricMatrix, TlVector>(DfObject::RUN_UKS_BETA);
        break;
    default:
        std::cerr << "program error. @DfConverge_Anderson::convergePMatrix()" << std::endl;
        break;
    }
}


