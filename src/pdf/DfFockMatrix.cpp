#include "DfFockMatrix.h"
#include "Fl_Geometry.h"
#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"
#include "TlVector.h"
#include "DfEriX.h"
#include "DfOverlapX.h"

DfFockMatrix::DfFockMatrix(TlSerializeData* pPdfParam) : DfObject(pPdfParam)
{
    this->isUseNewEngine_ = (*pPdfParam)["new_engine"].getBoolean();
}


DfFockMatrix::~DfFockMatrix()
{
}

void DfFockMatrix::DfFockMatrixMain()
{
    if (this->m_bDiskUtilization == false) {
        // DIRECT SCHEME
        switch (this->m_nMethodType) {
        case METHOD_RKS:
            this->mainDIRECT_RKS();
            break;
        case METHOD_UKS:
            this->mainDIRECT_UKS();
            break;
        case METHOD_ROKS:
            this->mainDIRECT_ROKS(); // ROKS
            break;
        default:
            // unsupported
            CnErr.abort("unsupported SCF type. stop.");
        }
    } else {
//     // FILE scheme
//     switch (this->m_nMethodType) {
//     case METHOD_RKS:
//       this->mainFILE(RUN_RKS);  // RKS
//       break;
//     case METHOD_UKS:
//       this->mainFILE(RUN_UKS_ALPHA);  // UKS alpha spin
//       this->mainFILE(RUN_UKS_BETA);  // UKS beta spin
//       break;
//     case METHOD_ROKS:
//       this->mainFILE_ROKS(); // ROKS
//       break;
//     default:
//       // unsupported
//       std::cerr << "unsupported SCF type. stop." << std::endl;
//       CnErr.abort();
//     }
        CnErr.abort("do coding. sorry!");
    }
}

void DfFockMatrix::mainDIRECT_RKS()
{
    this->mainDIRECT_RKS<TlSymmetricMatrix>();
}

void DfFockMatrix::mainDIRECT_UKS()
{
    this->mainDIRECT_UKS<TlSymmetricMatrix>();
}

void DfFockMatrix::mainDIRECT_ROKS()
{
    this->mainDIRECT_ROKS<TlMatrix, TlSymmetricMatrix, TlVector, DfEriX, DfOverlapX>();
}

void DfFockMatrix::setXC_RI(const RUN_TYPE nRunType, TlSymmetricMatrix& F)
{
    this->setXC_RI<TlSymmetricMatrix, TlVector, DfOverlapX>(nRunType, F);
}

void DfFockMatrix::setXC_DIRECT(const RUN_TYPE nRunType, TlSymmetricMatrix& F)
{
    this->setXC_DIRECT<TlSymmetricMatrix>(nRunType, F);
}

void DfFockMatrix::setCoulomb(const METHOD_TYPE nMethodType, TlSymmetricMatrix& F)
{
    TlSymmetricMatrix J(this->m_nNumOfAOs);
    if (this->J_engine_ == J_ENGINE_RI_J) {
        this->setCoulomb<TlSymmetricMatrix, TlVector, DfEriX>(nMethodType, J);
        // if (this->isUseNewEngine_ == true) {
        //     this->logger(" use new engine\n");
        //     this->setCoulomb<TlSymmetricMatrix, TlVector, DfEriX>(nMethodType, J);
        // } else {
        //     this->setCoulomb<TlSymmetricMatrix, TlVector, DfEri>(nMethodType, J);
        // }
        F += J;
        
        // update method
        if (this->m_nIteration > 1) {
            const TlSymmetricMatrix prevJ = DfObject::getJMatrix<TlSymmetricMatrix>(this->m_nIteration -1);
            J += prevJ;
        }

        this->saveJMatrix(this->m_nIteration, J);
    } else {
        J = this->getJMatrix<TlSymmetricMatrix>(this->m_nIteration);

        // update method
        if (this->m_nIteration > 1) {
            const TlSymmetricMatrix prevJ = DfObject::getJMatrix<TlSymmetricMatrix>(this->m_nIteration -1);
            J -= prevJ;
        }

        F += J;
    }
}

TlSymmetricMatrix DfFockMatrix::getFpqMatrix(const RUN_TYPE nRunType, const int nIteration)
{
    return (DfObject::getFpqMatrix<TlSymmetricMatrix>(nRunType, nIteration));
}

TlVector DfFockMatrix::getRho(const RUN_TYPE nRunType, const int nIteration)
{
    return (DfObject::getRho<TlVector>(nRunType, nIteration));
}

TlVector DfFockMatrix::getMyu(const RUN_TYPE nRunType, const int nIteration)
{
    return (DfObject::getMyu<TlVector>(nRunType, nIteration));
}

// void DfFockMatrix::saveFpqMatrix(const RUN_TYPE nRunType, const TlSymmetricMatrix& F)
// {
//     DfObject::saveFpqMatrix<TlSymmetricMatrix>(nRunType, this->m_nIteration, F);
// }

