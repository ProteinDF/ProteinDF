#include "DfOverlapX_Parallel.h"
#include "DfTaskCtrl_Parallel.h"
#include "TlCommunicate.h"
#include "TlSparseSymmetricMatrix.h"

DfOverlapX_Parallel::DfOverlapX_Parallel(TlSerializeData* pPdfParam)
    : DfOverlapX(pPdfParam) {
}


DfOverlapX_Parallel::~DfOverlapX_Parallel()
{
}


void DfOverlapX_Parallel::logger(const std::string& str) const
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfOverlapX::logger(str);
    }
}


DfTaskCtrl* DfOverlapX_Parallel::getDfTaskCtrlObject() const
{
    DfTaskCtrl* pDfTaskCtrl = new DfTaskCtrl_Parallel(this->pPdfParam_);
    return pDfTaskCtrl;
}


void DfOverlapX_Parallel::finalize(TlMatrix* pMtx)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    rComm.allReduce_SUM(*pMtx);
}


void DfOverlapX_Parallel::finalize(TlSymmetricMatrix* pMtx)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    rComm.allReduce_SUM(*pMtx);
}


void DfOverlapX_Parallel::finalize(TlVector* pVct)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    rComm.allReduce_SUM(*pVct);
}


void DfOverlapX_Parallel::getSpqD(TlDistributeSymmetricMatrix* pSpq)
{
    assert(pSpq != NULL);
    const index_type numOfAOs = this->m_nNumOfAOs;
    pSpq->resize(numOfAOs);
    TlSparseSymmetricMatrix tmpSpq(numOfAOs);

    const TlOrbitalInfo orbitalInfo((*(this->pPdfParam_))["coordinates"],
                                    (*(this->pPdfParam_))["basis_sets"]);
    this->calcOverlap(orbitalInfo, &tmpSpq);

    this->loggerTime(" finalize");
    pSpq->mergeSparseMatrix(tmpSpq);
}


void DfOverlapX_Parallel::getSabD(TlDistributeSymmetricMatrix* pSab)
{
    assert(pSab != NULL);
    const index_type numOfAuxDens = this->m_nNumOfAux;
    pSab->resize(numOfAuxDens);
    TlSparseSymmetricMatrix tmpSab(numOfAuxDens);
    
    const TlOrbitalInfo_Density orbitalInfo_Density((*(this->pPdfParam_))["coordinates"],
                                                    (*(this->pPdfParam_))["basis_sets_j"]);
    this->calcOverlap(orbitalInfo_Density, &tmpSab);
    
    this->loggerTime(" finalize");
    pSab->mergeSparseMatrix(tmpSab);
}
