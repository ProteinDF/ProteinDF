#include "CnError.h"
#include "DfJMatrix.h"
#include "DfEriX.h"
#include "TlSymmetricMatrix.h"
#include "TlUtils.h"

DfJMatrix::DfJMatrix(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam) {
}


DfJMatrix::~DfJMatrix()
{
}


void DfJMatrix::buildJ()
{
    const index_type numOfAOs = this->m_nNumOfAOs;
    TlSymmetricMatrix J(numOfAOs);

    switch (this->J_engine_) {
    case J_ENGINE_RI_J:
        this->log_.info("use RI_J engine");
        this->getJ_RI(&J);
        break;
    case J_ENGINE_CD:
        this->log_.info("use CD engine");
        this->getJ_CD(&J);
        break;
    default:
        this->log_.info("use conventional engine");
        this->getJ_conventional(&J);
        break;
    }

    DfObject::saveJMatrix(this->m_nIteration, J);
}


void DfJMatrix::getJ_conventional(TlSymmetricMatrix* pJ)
{
    TlSymmetricMatrix P;
    if (this->isUpdateMethod_ == true) {
        P = this->getDiffDensityMatrix<TlSymmetricMatrix>(RUN_RKS, this->m_nIteration);
    } else {
        P = this->getPpqMatrix<TlSymmetricMatrix>(RUN_RKS, this->m_nIteration -1);
    }

    DfEriX dfEri(this->pPdfParam_);
    dfEri.getJpq(P, pJ);

    if (this->isUpdateMethod_ == true) {
        if (this->m_nIteration > 1) {
            const TlSymmetricMatrix prevJ = this->getJMatrix<TlSymmetricMatrix>(this->m_nIteration -1);
            *pJ += prevJ;
        }
    }
}


void DfJMatrix::getJ_RI(TlSymmetricMatrix* pJ)
{
    TlVector Rho;
    switch (this->m_nMethodType) {
    case METHOD_RKS:
        Rho = DfObject::getRho<TlVector>(RUN_RKS, this->m_nIteration);
        break;
    case METHOD_UKS:
        Rho = DfObject::getRho<TlVector>(RUN_UKS_ALPHA, this->m_nIteration);
        Rho += DfObject::getRho<TlVector>(RUN_UKS_BETA, this->m_nIteration);
        break;
    default:
        this->log_.critical("this method is unsopported. sorry.");
        CnErr.abort();
        break;
    }

    if (this->isUpdateMethod_ == true) {
        if (this->m_nIteration > 1) {
            TlVector prevRho;
            switch (this->m_nMethodType) {
            case METHOD_RKS:
                prevRho = DfObject::getRho<TlVector>(RUN_RKS, this->m_nIteration -1);
                break;
            case METHOD_UKS:
                prevRho = DfObject::getRho<TlVector>(RUN_UKS_ALPHA, this->m_nIteration -1);
                prevRho += DfObject::getRho<TlVector>(RUN_UKS_BETA, this->m_nIteration -1);
                break;
            default:
                this->log_.critical("this method is unsopported. sorry.");
                CnErr.abort();
                break;
            }

            Rho -= prevRho;
        }
    }

    DfEriX dfEri(this->pPdfParam_);
    dfEri.getdeltaHpqA(Rho, *pJ);

    if (this->isUpdateMethod_ == true) {
        if (this->m_nIteration > 1) {
            const TlSymmetricMatrix prevJ = DfObject::getJMatrix<TlSymmetricMatrix>(this->m_nIteration -1);
            *pJ += prevJ;
        }
    }
}


void DfJMatrix::getJ_CD(TlSymmetricMatrix* pJ)
{
    const index_type numOfAOs = this->m_nNumOfAOs;

    const TlSymmetricMatrix P = this->getPpqMatrix<TlSymmetricMatrix>(RUN_RKS, this->m_nIteration -1);

    // cholesky vector
    TlMatrix L;
    L.load("L.mat");
    const index_type numOfCBs = L.getNumOfCols();
    
    for (index_type I = 0; I < numOfCBs; ++I) {
        TlSymmetricMatrix LI(numOfAOs);
        for (index_type p = 0; p < numOfAOs; ++p) {
            for (index_type q = 0; q <= p; ++q) {
                const index_type rs = p + (2 * numOfAOs - (q +1)) * q / 2;
                LI.set(p, q, L.get(rs, I));
            }
        }

        TlMatrix QI = LI;
        QI.dot(P);
        double qi = QI.sum();

        *pJ += qi*LI;
    }
}
