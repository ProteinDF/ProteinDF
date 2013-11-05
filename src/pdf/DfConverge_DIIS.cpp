#include <cassert>
#include "DfConverge_DIIS.h"
#include "TlSymmetricMatrix.h"
#include "TlTime.h"

DfConverge_DIIS::DfConverge_DIIS(TlSerializeData* pPdfParam)
    : DfConverge_Damping(pPdfParam)
{
    const TlSerializeData& pdfParam = *pPdfParam;
    this->startIterationOfDIIS_ =
        std::max(pdfParam["scf_acceleration/DIIS/start"].getInt(), 1);
    this->numOfLastItems_ =
        std::max(pdfParam["scf_acceleration/DIIS/last_items"].getInt(), 3);

    if (this->pPdfParam_->hasKey("DIIS_e_max_pass_itr") == false) {
        (*this->pPdfParam_)["DIIS_e_max_pass_itr"] = -1;
    }
}

DfConverge_DIIS::~DfConverge_DIIS()
{
}

void DfConverge_DIIS::convergeRhoTilde()
{
    this->log_.critical("NOT support: DIIS for rho.");

    // switch (this->m_nMethodType) {
    // case METHOD_RKS:
    //     this->convergeRhoTilde<TlMatrix, TlSymmetricMatrix, TlVector>(DfObject::RUN_RKS);
    //     break;
    // case METHOD_UKS:
    //     this->convergeRhoTilde<TlMatrix, TlSymmetricMatrix, TlVector>(DfObject::RUN_UKS_ALPHA);
    //     this->convergeRhoTilde<TlMatrix, TlSymmetricMatrix, TlVector>(DfObject::RUN_UKS_BETA);
    //     break;
    // case METHOD_ROKS:
    //     this->convergeRhoTilde<TlMatrix, TlSymmetricMatrix, TlVector>(DfObject::RUN_UKS_ALPHA);
    //     this->convergeRhoTilde<TlMatrix, TlSymmetricMatrix, TlVector>(DfObject::RUN_UKS_BETA);
    //     break;
    // default:
    //     std::cerr << "program error. @DfConverge_Damping::convergeRhoTilde()" << std::endl;
    //     break;
    // }
}

void DfConverge_DIIS::convergeKSMatrix()
{
    switch (this->m_nMethodType) {
    case METHOD_RKS:
        this->convergeKSMatrix<TlMatrix, TlSymmetricMatrix>(DfObject::RUN_RKS);
        break;
    case METHOD_UKS:
        this->convergeKSMatrix<TlMatrix, TlSymmetricMatrix>(DfObject::RUN_UKS_ALPHA);
        this->convergeKSMatrix<TlMatrix, TlSymmetricMatrix>(DfObject::RUN_UKS_BETA);
        break;
    case METHOD_ROKS:
        this->convergeKSMatrix<TlMatrix, TlSymmetricMatrix>(DfObject::RUN_UKS_ALPHA);
        this->convergeKSMatrix<TlMatrix, TlSymmetricMatrix>(DfObject::RUN_UKS_BETA);
        break;
    default:
        std::cerr << "program error. @DfConverge_Damping::convergeKSMatrix()" << std::endl;
        break;
    }
}

void DfConverge_DIIS::convergePMatrix()
{
    switch (this->m_nMethodType) {
    case METHOD_RKS:
        this->convergePMatrix<TlMatrix, TlSymmetricMatrix>(DfObject::RUN_RKS);
        break;
    case METHOD_UKS:
        this->convergePMatrix<TlMatrix, TlSymmetricMatrix>(DfObject::RUN_UKS_ALPHA);
        this->convergePMatrix<TlMatrix, TlSymmetricMatrix>(DfObject::RUN_UKS_BETA);
        break;
    case METHOD_ROKS:
        this->convergePMatrix<TlMatrix, TlSymmetricMatrix>(DfObject::RUN_UKS_ALPHA);
        this->convergePMatrix<TlMatrix, TlSymmetricMatrix>(DfObject::RUN_UKS_BETA);
        break;
    default:
        std::cerr << "program error. @DfConverge_Damping::convergePMatrix()" << std::endl;
        break;
    }
}

