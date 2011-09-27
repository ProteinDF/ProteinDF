#include "DfOverlap_Parallel.h"

#include "TlFile.h"
#include "TlCommunicate.h"
#include "TlSparseSymmetricMatrix.h"
#include "TlFileSymmetricMatrix.h"

DfOverlap_Parallel::DfOverlap_Parallel(TlSerializeData* pPdfParam)
    : DfOverlap(pPdfParam)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    rComm.barrier();
}


DfOverlap_Parallel::~DfOverlap_Parallel()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    rComm.barrier();
}


void DfOverlap_Parallel::logger(const std::string& str) const
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    if (rComm.isMaster() == true) {
        DfOverlap::logger(str);
    }
}


void DfOverlap_Parallel::cutoffReport()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    rComm.allReduce_SUM(this->CcountSS);
    rComm.allReduce_SUM(this->CcountPS);
    rComm.allReduce_SUM(this->CcountPP);
    rComm.allReduce_SUM(this->CcountDS);
    rComm.allReduce_SUM(this->CcountDP);
    rComm.allReduce_SUM(this->CcountDD);
    rComm.allReduce_SUM(this->TcountSS);
    rComm.allReduce_SUM(this->TcountPS);
    rComm.allReduce_SUM(this->TcountPP);
    rComm.allReduce_SUM(this->TcountDS);
    rComm.allReduce_SUM(this->TcountDP);
    rComm.allReduce_SUM(this->TcountDD);

    DfOverlap::cutoffReport();
}


std::vector<DfObject::IJShellPair>
DfOverlap_Parallel::getQueue(const TlOrbitalInfo* pOrbitalInfo,
                             const std::vector<std::vector<std::size_t> >& shellList,
                             const int maxScore, const bool initialize)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int nProc = rComm.getNumOfProc();
    const int nRank = rComm.getRank();

    static std::size_t count = 0;
    if (initialize == true) {
        count = 0;
    }

    std::vector<IJShellPair> ij = DfOverlap::getQueue(pOrbitalInfo, shellList,
                                                      maxScore, initialize);
    std::vector<IJShellPair> ij_part;
    for (std::vector<IJShellPair>::const_iterator p = ij.begin(); p != ij.end(); ++p) {
        if ((count % nProc) == static_cast<std::size_t>(nRank)) {
            ij_part.push_back(*p);
        }
        ++count;
    }
    
    return ij_part;
}


std::vector<DfObject::IJShellPair>
DfOverlap_Parallel::getQueue(const TlOrbitalInfo& orbitalInfo1,
                             const TlOrbitalInfo& orbitalInfo2,
                             const std::vector<std::vector<std::size_t> >& shellList1,
                             const std::vector<std::vector<std::size_t> >& shellList2,
                             const int maxNumOfPairs, const bool initialize)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int nProc = rComm.getNumOfProc();
    const int nRank = rComm.getRank();

    static std::size_t count = 0;
    if (initialize == true) {
        count = 0;
    }

    std::vector<IJShellPair> ij = DfOverlap::getQueue(orbitalInfo1,
                                                      orbitalInfo2,
                                                      shellList1, shellList2,
                                                      maxNumOfPairs, initialize);
    std::vector<IJShellPair> ij_part;
    for (std::vector<IJShellPair>::const_iterator p = ij.begin(); p != ij.end(); ++p) {
        if ((count % nProc) == static_cast<std::size_t>(nRank)) {
            ij_part.push_back(*p);
        }
        ++count;
    }

    return ij_part;
}


// 積分の後処理
void DfOverlap_Parallel::finalizeIntegral(TlSymmetricMatrix& rMatrix)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    rComm.allReduce_SUM(rMatrix);
}


void DfOverlap_Parallel::finalizeIntegral(TlVector& rVector)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    rComm.allReduce_SUM(rVector);
}

// =====================================================================
// for LAPACK
// void DfOverlap_Parallel::getSgd(TlSymmetricMatrix& Sgd) {
//   DfOverlap::getSgd(Sgd);
// }

void DfOverlap_Parallel::getdeltaHpqG(const TlVector& deltaMyu, TlSymmetricMatrix& deltaH)
{
    DfOverlap::getdeltaHpqG(deltaMyu, deltaH);
}

// =====================================================================
// for ScaLAPACK
// void DfOverlap_Parallel::getSgd(TlDistributeSymmetricMatrix& Sgd) {
//   this->loggerTime(" start");

//   this->loggerTime(" make table");
//   this->makeTable();

//   //this->initializeIntegral();
//   TlSparseSymmetricMatrix SgdTmp(Sgd.getNumOfRows());

//   this->loggerTime(" integral");
//   const int maxScore = this->blockSize_;
//   std::list<IJShellPair> ij = this->getQueueForAux(maxScore, true);
//   while (ij.empty() != true) {
//     this->sgdcalc(ij, &SgdTmp);
//     ij = this->getQueueForAux(maxScore);
//   }

//   this->loggerTime(" end");
//   Sgd.mergePartialMatrix(SgdTmp);
// }


void DfOverlap_Parallel::getdeltaHpqG(const TlVector& deltaMyu,
                                      TlDistributeSymmetricMatrix& deltaH)
{
    this->loggerTime(" start");

    // Set new cutvalue
    const double MAXdeltaMyu = deltaMyu.getMaxAbsoluteElement();
    if ((MAXdeltaMyu > cutvalue) && (MAXdeltaMyu < 1.0)) {
        cutvalue /= MAXdeltaMyu;
        this->logger(TlUtils::format("new cutvalue in DfOverlap.getdeltaHpqG is %.2e\n", cutvalue));
    }

    this->loggerTime(" make table");
    this->makeTable();

    TcountSS = TcountPS = TcountPP = TcountDS = TcountDP = TcountDD = 0;
    CcountSS = CcountPS = CcountPP = CcountDS = CcountDP = CcountDD = 0;

    //this->initializeIntegral();
    TlSparseSymmetricMatrix deltaHTmp(deltaH.getNumOfRows());

    this->loggerTime(" integral");
    const int maxScore = this->blockSize_;
    std::vector<IJShellPair> ij = this->getQueue(this->pOrbitalInfo_,
                                                 this->shellList_, maxScore, true);
    while (ij.empty() != true) {
        this->ovpqgDH(ij, deltaMyu, &deltaHTmp);
        ij = this->getQueue(this->pOrbitalInfo_, this->shellList_, maxScore);
    }

    this->loggerTime("finalize");
    this->cutoffReport();
    deltaH.mergeSparseMatrix(deltaHTmp);
    this->loggerTime(" end");
}


void DfOverlap_Parallel::getdeltaHpqG(const TlDistributeVector& deltaMyu,
                                      TlDistributeSymmetricMatrix& deltaH)
{
    this->loggerTime(" gather distribute vector.");
    const TlVector globalDeltaMyu(deltaMyu.getVector());
    this->getdeltaHpqG(globalDeltaMyu, deltaH);
}


// [alpha]
// <p q> ===============================================================
void DfOverlap_Parallel::getNa(TlVector* pNa)
{
    assert(pNa != NULL);
    pNa->resize(this->m_nNumOfAux);

    this->getNa_core(pNa);
}


// <p q> ===============================================================
void DfOverlap_Parallel::getSpq(TlSymmetricMatrix* pSpq)
{
    pSpq->resize(this->m_nNumOfAOs);

    this->getSpq_core(pSpq);

    this->loggerTime(" finalize");
    this->cutoffReport();

    TlCommunicate& rComm = TlCommunicate::getInstance();
    rComm.allReduce_SUM(*pSpq);
}


void DfOverlap_Parallel::getSpq(TlDistributeSymmetricMatrix* pSpq)
{
    pSpq->resize(this->m_nNumOfAOs);

    TlSparseSymmetricMatrix Spq_part(this->m_nNumOfAOs);
    this->getSpq_core(&Spq_part);

    this->loggerTime(" finalize");
    this->cutoffReport();
    pSpq->mergeSparseMatrix(Spq_part);
}


// <alpha beta> ========================================================
void DfOverlap_Parallel::getSab2(TlSymmetricMatrix* pSab2)
{
    assert(pSab2 != NULL);
    pSab2->resize(this->m_nNumOfAux);
    
    this->getSab2_core(pSab2);

    this->loggerTime(" finalize");
    TlCommunicate& rComm = TlCommunicate::getInstance();
    rComm.allReduce_SUM(*pSab2);
}


void DfOverlap_Parallel::getSab2(TlDistributeSymmetricMatrix* pSab2)
{
    assert(pSab2 != NULL);
    pSab2->resize(this->m_nNumOfAux);

    TlSparseSymmetricMatrix Sab2_part(this->m_nNumOfAux);
    this->getSab2_core(&Sab2_part);

    this->loggerTime(" finalize");
    pSab2->mergeSparseMatrix(Sab2_part);
}


// <gamma delta> ========================================================
void DfOverlap_Parallel::getSgd(TlSymmetricMatrix* pSgd)
{
    assert(pSgd != NULL);
    pSgd->resize(this->numOfAuxXC_);
    
    this->getSgd_core(pSgd);

    this->loggerTime(" finalize");
    TlCommunicate& rComm = TlCommunicate::getInstance();
    rComm.allReduce_SUM(*pSgd);
}


void DfOverlap_Parallel::getSgd(TlDistributeSymmetricMatrix* pSgd)
{
    assert(pSgd != NULL);
    pSgd->resize(this->numOfAuxXC_);

    TlSparseSymmetricMatrix Sgd_part(this->numOfAuxXC_);
    this->getSgd_core(&Sgd_part);

    this->loggerTime(" finalize");
    pSgd->mergeSparseMatrix(Sgd_part);
}


