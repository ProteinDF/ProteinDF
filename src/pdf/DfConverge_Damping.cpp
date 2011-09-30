#include <cassert>
#include "DfConverge_Damping.h"
#include "TlSymmetricMatrix.h"
#include "TlLogX.h"
#include "TlTime.h"

DfConverge_Damping::DfConverge_Damping(TlSerializeData* pPdfParam)
    : DfConverge(pPdfParam)
{
    const TlSerializeData& pdfParam = *pPdfParam;
    this->m_nStartIteration =
        std::max(pdfParam["scf-acceleration/damping/number-of-damping"].getInt(), 2);
    this->m_dDampingFactor =
        pdfParam["scf-acceleration/damping/damping-factor"].getDouble();
}

DfConverge_Damping::~DfConverge_Damping()
{
}

void DfConverge_Damping::convergeRhoTilde()
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
        std::cerr << "program error. @DfConverge_Damping::convergeRhoTilde()" << std::endl;
        break;
    }
}

void DfConverge_Damping::convergeKSMatrix()
{
    switch (this->m_nMethodType) {
    case METHOD_RKS:
        this->convergeKSMatrix<TlSymmetricMatrix>(DfObject::RUN_RKS);
        break;
    case METHOD_UKS:
        this->convergeKSMatrix<TlSymmetricMatrix>(DfObject::RUN_UKS_ALPHA);
        this->convergeKSMatrix<TlSymmetricMatrix>(DfObject::RUN_UKS_BETA);
        break;
    case METHOD_ROKS:
        this->convergeKSMatrix<TlSymmetricMatrix>(DfObject::RUN_UKS_ALPHA);
        this->convergeKSMatrix<TlSymmetricMatrix>(DfObject::RUN_UKS_BETA);
        break;
    default:
        std::cerr << "program error. @DfConverge_Damping::convergeKSMatrix()" << std::endl;
        break;
    }
}

void DfConverge_Damping::convergePMatrix()
{
    switch (this->m_nMethodType) {
    case METHOD_RKS:
        this->convergePMatrix<TlSymmetricMatrix>(DfObject::RUN_RKS);
        break;
    case METHOD_UKS:
        this->convergePMatrix<TlSymmetricMatrix>(DfObject::RUN_UKS_ALPHA);
        this->convergePMatrix<TlSymmetricMatrix>(DfObject::RUN_UKS_BETA);
        break;
    case METHOD_ROKS:
        this->convergePMatrix<TlSymmetricMatrix>(DfObject::RUN_UKS_ALPHA);
        this->convergePMatrix<TlSymmetricMatrix>(DfObject::RUN_UKS_BETA);
        break;
    default:
        std::cerr << "program error. @DfConverge_Damping::convergePMatrix()" << std::endl;
        break;
    }
}

