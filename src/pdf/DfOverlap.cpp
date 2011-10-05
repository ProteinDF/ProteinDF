#include <cstdlib>
#include <cmath>
#include <cassert>

#include "DfOverlap.h"
#include "TlTime.h"
#include "TlMatrix.h"
#include "TlLogX.h"

#define MAX_SHELL 2 // s=0, p=1, d=2
#define DEFAULT_BLOCK_SIZE 1000

const double DfOverlap::SQR3I  = 1.0 / sqrt(3.0);
const double DfOverlap::SQR3I2 = 2.0 / sqrt(3.0);

DfOverlap::DfOverlap(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam), pOrbitalInfo_(NULL),
      pOrbitalInfo_Density_(NULL), pOrbitalInfo_XC_(NULL)
{
    this->pOrbitalInfo_ = new TlOrbitalInfo((*pPdfParam)["coordinates"],
                                            (*pPdfParam)["basis_sets"]);
    this->pOrbitalInfo_Density_ = new TlOrbitalInfo_Density((*pPdfParam)["coordinates"],
                                                            (*pPdfParam)["basis_sets_j"]);
    this->pOrbitalInfo_XC_ = new TlOrbitalInfo_XC((*pPdfParam)["coordinates"],
                                                  (*pPdfParam)["basis_sets_k"]);

    const TlSerializeData& pdfParam = *pPdfParam;
    this->cutvalue = pdfParam["cut-value"].getDouble();

    this->blockSize_ = DEFAULT_BLOCK_SIZE;
    if (pdfParam["block-size"].getStr().empty() != true) {
        this->blockSize_ = pdfParam["block-size"].getInt();
    }
}



DfOverlap::~DfOverlap()
{
    if (this->pOrbitalInfo_ != NULL) {
        delete this->pOrbitalInfo_;
    }

    if (this->pOrbitalInfo_Density_ != NULL) {
        delete this->pOrbitalInfo_Density_;
    }

    if (this->pOrbitalInfo_XC_ != NULL) {
        delete this->pOrbitalInfo_XC_;
    }
}


// 統計用
void DfOverlap::countupTotal(int ity, int jty)
{
    const int nType = ity * 4 + jty;
    // SS = 0,
    // PS = 4, PP = 5,
    // DS = 8, DP = 9, DD = 10
    switch (nType) {
    case 0:
        ++(this->TcountSS);
        break;
    case 4:
        ++(this->TcountPS);
        break;
    case 5:
        ++(this->TcountPP);
        break;
    case 8:
        ++(this->TcountDS);
        break;
    case 9:
        ++(this->TcountDP);
        break;
    case 10:
        ++(this->TcountDD);
        break;
    }
}

// cutoffした回数をカウント
void DfOverlap::countupCutoff(int ity, int jty)
{
    const int nType = ity * 4 + jty;
    switch (nType) {
    case 0:
        ++(this->CcountSS);
        break;
    case 4:
        ++(this->CcountPS);
        break;
    case 5:
        ++(this->CcountPP);
        break;
    case 8:
        ++(this->CcountDS);
        break;
    case 9:
        ++(this->CcountDP);
        break;
    case 10:
        ++(this->CcountDD);
        break;
    }
}


void DfOverlap::cutoffReport()
{
    this->logger(" cut off report:\n");
    this->cutoffReport("ss", this->CcountSS, this->TcountSS);
    this->cutoffReport("ps", this->CcountPS, this->TcountPS);
    this->cutoffReport("pp", this->CcountPP, this->TcountPP);
    this->cutoffReport("ds", this->CcountDS, this->TcountDS);
    this->cutoffReport("dp", this->CcountDP, this->TcountDP);
    this->cutoffReport("dd", this->CcountDD, this->TcountDD);
    this->logger("\n");
}


void DfOverlap::cutoffReport(const std::string& shell,
                             const std::size_t cutoffCount, const std::size_t totalCount)
{
    double rate = 0.0;
    if (totalCount > 0) {
        rate = (double)cutoffCount / (double)totalCount * 100.0;
    }

    this->logger(TlUtils::format(" shell %s: %12lu/%12lu (%6.2f%%)\n",
                                 shell.c_str(), cutoffCount, totalCount, rate));
}


void DfOverlap::getNa(TlVector* pNa)
{
    assert(pNa != NULL);
    pNa->resize(this->m_nNumOfAux);

    this->getNa_core(pNa);
}


void DfOverlap::getNa_core(TlVector* pNa)
{
    assert(pNa != NULL);

    this->loggerTime(" make aux table");
    this->makeTable();

    this->loggerTime(" integral");

    const int MAXTYPE = 3;
    for (int kty = 0; kty < MAXTYPE; ++kty) {
//         const std::size_t max_k = this->shellList_Dens_[kty].size();
//         for (std::size_t k = 0; k < max_k; ++k) {
//             calcList.push_back(this->shellList_Dens_[kty][k]);
//         }
        const std::vector<std::size_t> calcList = this->shellList_Dens_[kty];

        this->calcNa(calcList, pNa);
    }
}


// [alpha beta]
void DfOverlap::getSab2(TlSymmetricMatrix* pSab2)
{
    assert(pSab2 != NULL);

    pSab2->resize(this->m_nNumOfAux);
    this->getSab2_core(pSab2);

    this->loggerTime(" finalize");
    //this->cutoffReport();
}


void DfOverlap::getSab2_core(TlMatrixObject* pSab2)
{
    assert(pSab2 != NULL);

    this->loggerTime(" make aux table");
    this->makeTable();

    this->loggerTime(" integral");
    const std::size_t maxNumOfPairs = this->blockSize_;
    std::vector<IJShellPair> jk = this->getQueue(NULL, this->shellList_Dens_,
                                                 maxNumOfPairs, true);
    while (jk.empty() != true) {
        this->spqcalc(*(this->pOrbitalInfo_Density_), jk, pSab2);
        jk = this->getQueue(NULL, this->shellList_Dens_, maxNumOfPairs);
    }
}


// [gamma delta]
void DfOverlap::getSgd(TlSymmetricMatrix* pSgd)
{
    assert(pSgd != NULL);
    pSgd->resize(this->numOfAuxXC_);

    this->getSgd_core(pSgd);
}


void DfOverlap::getSgd_core(TlMatrixObject* pSgd)
{
    assert(pSgd != NULL);

    this->loggerTime(" make table");
    this->makeTable();

    this->loggerTime(" integral");
    const int maxNumOfPairs = this->blockSize_ / sizeof(double);
    std::vector<IJShellPair> jk = this->getQueue(NULL, this->shellList_XC_, maxNumOfPairs, true);
    while (jk.empty() != true) {
        this->spqcalc(*(this->pOrbitalInfo_XC_), jk, pSgd);
        jk = this->getQueue(NULL, this->shellList_XC_, maxNumOfPairs);
    }
}


////////////////////////////////////////////////////////////////////////
void DfOverlap::getSpq(TlSymmetricMatrix* pSpq)
{
    assert(pSpq != NULL);

    pSpq->resize(this->m_nNumOfAOs);
    this->getSpq_core(pSpq);

    this->loggerTime(" finalize");
    this->cutoffReport();
}


void DfOverlap::getSpq_core(TlMatrixObject* pSpq)
{
    assert(pSpq != NULL);

    this->loggerTime(" make table");
    this->makeTable();

    // initialize counters
    TcountSS = TcountPS = TcountPP = TcountDS = TcountDP = TcountDD = 0;
    CcountSS = CcountPS = CcountPP = CcountDS = CcountDP = CcountDD = 0;

    this->loggerTime(" integral");
    const std::size_t maxNumOfPairs = this->blockSize_ / sizeof(double);
    std::vector<IJShellPair> ij = this->getQueue(this->pOrbitalInfo_,
                                                 this->shellList_,
                                                 maxNumOfPairs, true);
    while (ij.empty() != true) {
        this->spqcalc(*(this->pOrbitalInfo_), ij, pSpq);
        ij = this->getQueue(this->pOrbitalInfo_, this->shellList_, maxNumOfPairs);
    }
}


// void DfOverlap::getSpq_new()
// {
//     this->loggerTime(" start (new Overlap)");
//     TlSymmetricMatrix Spq(this->m_nNumOfAOs);

//     this->getSpq_core2(&Spq);

//     this->loggerTime(" finalize");
//     Spq.save("newSpq.matrix");
// }


// void DfOverlap::getSpq_core2(TlMatrixObject* pSpq)
// {
//     assert(pSpq != NULL);

//     // initialize counters
//     TcountSS = TcountPS = TcountPP = TcountDS = TcountDP = TcountDD = 0;
//     CcountSS = CcountPS = CcountPP = CcountDS = CcountDP = CcountDD = 0;
//      DfOvp ovp;
    
//     this->loggerTime(" integral");
//     const std::size_t maxNumOfPairs = this->blockSize_ / sizeof(double);
//     std::vector<IJShellPair> ij = this->getQueue(this->pOrbitalInfo_,
//                                                  this->shellList_,
//                                                  maxNumOfPairs, true);
//     const TlPosition C(0.0, 0.0, 0.0);
//     DfOvp::PGTOs pGTOs_C(1);
//     pGTOs_C[0] = DfOvp::PGTO(1.0, 0.0);
//     while (ij.empty() != true) {
//         for (std::vector<IJShellPair>::const_iterator p = ij.begin(); p != ij.end(); ++p) {
//             const int shellIndexA = p->nIShell;
//             const int shellIndexB = p->nJShell;
//             const int shellTypeA = this->pOrbitalInfo_->getShellType(shellIndexA);
//             const int shellTypeB = this->pOrbitalInfo_->getShellType(shellIndexB);
//             const DfOvp::Query query(0, 0, 0, shellTypeA, shellTypeB, 0);
//             const TlPosition A = this->pOrbitalInfo_->getPosition(shellIndexA);
//             const TlPosition B = this->pOrbitalInfo_->getPosition(shellIndexB);

//             const int numOfPGTOsA = this->pOrbitalInfo_->getCgtoContraction(shellIndexA);
//             DfOvp::PGTOs pGTOs_A(numOfPGTOsA);
//             for (int i = 0; i < numOfPGTOsA; ++i) {
//                 const DfOvp::PGTO pgto(this->pOrbitalInfo_->getCoefficient(shellIndexA, i),
//                                        this->pOrbitalInfo_->getExponent(shellIndexA, i));
//                 pGTOs_A[i] = pgto;
//             }

//             const int numOfPGTOsB = this->pOrbitalInfo_->getCgtoContraction(shellIndexB);
//             DfOvp::PGTOs pGTOs_B(numOfPGTOsB);
//             for (int i = 0; i < numOfPGTOsB; ++i) {
//                 const DfOvp::PGTO pgto(this->pOrbitalInfo_->getCoefficient(shellIndexB, i),
//                                        this->pOrbitalInfo_->getExponent(shellIndexB, i));
//                 pGTOs_B[i] = pgto;
//             }

//             ovp.calc(query, A, B, C, pGTOs_A, pGTOs_B, pGTOs_C);

//             const int dimA = shellTypeA * 2 + 1;
//             const int dimB = shellTypeB * 2 + 1;
//             int index = 0;
//             for (int i = 0; i < dimA; ++i) {
//                 const index_type I = shellIndexA + i;
//                 for (int j = 0; j < dimB; ++j) {
//                     const index_type J = shellIndexB + j;

//                     const double value = ovp.WORK[index];
//                     pSpq->set(I, J, value);
//                     ++index;
//                 }
//             }
//             std::cerr << std::endl;
//         }
        
//         ij = this->getQueue(this->pOrbitalInfo_, this->shellList_, maxNumOfPairs);
//     }
// }


TlMatrix DfOverlap::getSpq(const TlOrbitalInfo& orbInfo1, const TlOrbitalInfo& orbInfo2)
{
    //this->loggerTime(" start");
    const std::size_t numOfAOs1 = orbInfo1.getNumOfOrbitals();
    const std::size_t numOfAOs2 = orbInfo2.getNumOfOrbitals();
    TlMatrix Spq(numOfAOs1, numOfAOs2);

    this->getSpq_core(orbInfo1, orbInfo2, &Spq);

    //this->loggerTime(" finalize");
    //this->cutoffReport();
    //Spq.save(savePath);
    
    //this->loggerTime(" end");
    return Spq;
}


void DfOverlap::getSpq_core(const TlOrbitalInfo& orbInfo1,
                            const TlOrbitalInfo& orbInfo2,
                            TlMatrixObject* pSpq)
{
    assert(pSpq != NULL);

    //this->loggerTime(" make table");
    const ShellListType shellList1 = this->makeShellList(orbInfo1);
    const ShellListType shellList2 = this->makeShellList(orbInfo2);

    // initialize counters
    TcountSS = TcountPS = TcountPP = TcountDS = TcountDP = TcountDD = 0;
    CcountSS = CcountPS = CcountPP = CcountDS = CcountDP = CcountDD = 0;

    this->loggerTime(" integral");
    const std::size_t maxNumOfPairs = this->blockSize_ / sizeof(double);
    std::vector<IJShellPair> ij = this->getQueue(orbInfo1,
                                                 orbInfo2,
                                                 shellList1,
                                                 shellList2,
                                                 maxNumOfPairs, true);
    while (ij.empty() != true) {
        this->spqcalc(orbInfo1, orbInfo2, ij, pSpq);
        ij = this->getQueue(orbInfo1,
                            orbInfo2,
                            shellList1,
                            shellList2,
                            maxNumOfPairs);
    }
}


////////////////////////////////////////////////////////////////////////
void DfOverlap::getdeltaHpqG(const TlVector& deltaMyu, TlSymmetricMatrix& deltaH)
{
    this->loggerTime(" start");

    // Set new cutvalue
    const double MAXdeltaMyu = deltaMyu.getMaxAbsoluteElement();
    if ((MAXdeltaMyu > cutvalue) && (MAXdeltaMyu < 1.0)) {
        cutvalue /= MAXdeltaMyu;
        this->logger(TlUtils::format("new cutvalue is %.2e\n", cutvalue));
    }

    this->loggerTime(" make table");
    this->makeTable();

    TcountSS = TcountPS = TcountPP = TcountDS = TcountDP = TcountDD = 0;
    CcountSS = CcountPS = CcountPP = CcountDS = CcountDP = CcountDD = 0;

    //this->initializeIntegral();

    this->loggerTime(" integral");
    const std::size_t maxNumOfPairs = this->blockSize_;
    std::vector<IJShellPair> ij = this->getQueue(this->pOrbitalInfo_,
                                                 this->shellList_,
                                                 maxNumOfPairs, true);
    while (ij.empty() != true) {
        this->ovpqgDH(ij, deltaMyu, &deltaH);
        ij = this->getQueue(this->pOrbitalInfo_, this->shellList_, maxNumOfPairs);
    }

    this->loggerTime(" finalize");
    this->finalizeIntegral(deltaH);
    this->cutoffReport();

    this->loggerTime(" end");
}

void DfOverlap::getdeltaHpqG(const TlVector& deltaMyu, const TlVector& deltaEps,
                             TlSymmetricMatrix& deltaH1, TlSymmetricMatrix& deltaH2)
{
    this->loggerTime(" start");

    const double MAXdeltaMyu = deltaMyu.getMaxAbsoluteElement();
    const double MAXdeltaEps = deltaEps.getMaxAbsoluteElement();

    double cutvalue_F = this->cutvalue;
    double cutvalue_E = this->cutvalue;
    if ((MAXdeltaMyu > cutvalue_F) && (MAXdeltaMyu < 1.0)) {
        cutvalue_F /= MAXdeltaMyu;
    }
    if ((MAXdeltaEps > cutvalue_E) && (MAXdeltaEps < 1.0)) {
        cutvalue_E /= MAXdeltaEps;
    }

    this->cutvalue = std::max(cutvalue_F, cutvalue_E);
//   log << "new cutvalue in DfOverlap.getdeltaHpqG is "   << cutvalue << "\n";
//   log << "new cutvalue_F in DfOverlap.getdeltaHpqG is " << cutvalue_F << "\n";
//   log << "new cutvalue_E in DfOverlap.getdeltaHpqG is " << cutvalue_E << "\n";

    this->loggerTime(" make table");
    this->makeTable();

    TcountSS = TcountPS = TcountPP = TcountDS = TcountDP = TcountDD = 0;
    CcountSS = CcountPS = CcountPP = CcountDS = CcountDP = CcountDD = 0;

    //this->initializeIntegral();

    this->loggerTime(" integral");
    const std::size_t maxNumOfPairs = this->blockSize_;
    std::vector<IJShellPair> ij = this->getQueue(this->pOrbitalInfo_,
                                                 this->shellList_, maxNumOfPairs, true);
    while (ij.empty() != true) {
        this->ovpqgDH(ij, deltaMyu, deltaEps, &deltaH1, &deltaH2);
        ij = this->getQueue(this->pOrbitalInfo_, this->shellList_, maxNumOfPairs);
    }

    this->loggerTime(" finalize");
    this->finalizeIntegral(deltaH1);
    this->finalizeIntegral(deltaH2);

    this->cutoffReport();

    this->loggerTime(" end");
}


std::vector<DfObject::IJShellPair>
DfOverlap::getQueue(const TlOrbitalInfo* pOrbitalInfo,
                    const std::vector<std::vector<std::size_t> >& shellList,
                    const int maxNumOfPairs, const bool initialize)
{
    static const int MAXTYPE = 3; // s, p, d
    static bool hasCalced = false;
    struct State {
        int ity;
        int jty;
        int i;
        int j;
    };
    static State state;
    if (initialize == true) {
        hasCalced = false;
        state.ity = MAXTYPE -1;
        state.jty = state.ity;
        state.i = -1;
        state.j = -1;
    }

    std::vector<IJShellPair> ijShellPairList;
    ijShellPairList.reserve(maxNumOfPairs);
    if (hasCalced == true) {
        return ijShellPairList;
    }

    int numOfPairs = 0;
    // I-shell Type ------------------------------------------------------
    for (int ity = state.ity; ity >= 0; --ity) {

        // J-shell Type ----------------------------------------------------
        for (int jty = std::min(ity, state.jty); jty >= 0; --jty) {

            for (std::size_t i = std::max(0, state.i); i < shellList[ity].size(); ++i) {
                const std::size_t index_i = shellList[ity][i];

                for (std::size_t j = std::max(0, state.j); j < shellList[jty].size(); ++j) {
                    const std::size_t index_j = shellList[jty][j];

                    if ((ity == jty) && (index_i < index_j)) {
                        continue;
                    }

                    // cut off
                    if (pOrbitalInfo != NULL) {
                        this->countupTotal(ity, jty);
                        if (this->isCutoffUsingSchwartzInequality(*pOrbitalInfo, index_i, index_j) == true) {
                            this->countupCutoff(ity, jty);
                            continue;
                        }
                    }

                    ijShellPairList.push_back(IJShellPair(index_i,
                                                          index_j));
                    //score += (ity +1) * (jty +1);
                    ++numOfPairs;

                    if (numOfPairs > maxNumOfPairs) {
                        // 現状の保存
                        state.ity = ity;
                        state.jty = jty;
                        state.i = i;
                        state.j = j +1; // "+1" means to start next j-shell at next time.

                        return ijShellPairList;
                    }
                }
                state.j = -1;
            }
            state.i = -1;
        }
        state.jty = MAXTYPE;
    }

    // finalize
    hasCalced = true;
    return ijShellPairList;
}

std::vector<DfObject::IJShellPair>
DfOverlap::getQueue(const TlOrbitalInfo& orbitalInfo1,
                    const TlOrbitalInfo& orbitalInfo2,
                    const std::vector<std::vector<std::size_t> >& shellList1,
                    const std::vector<std::vector<std::size_t> >& shellList2,
                    const int maxNumOfPairs, const bool initialize)
{
    static const int MAXTYPE = 3; // s, p, d
    static bool hasCalced = false;
    struct State {
        int ity;
        int jty;
        int i;
        int j;
    };
    static State state;
    if (initialize == true) {
        hasCalced = false;
        state.ity = MAXTYPE -1;
        state.jty = state.ity;
        state.i = -1;
        state.j = -1;
    }

    std::vector<IJShellPair> ijShellPairList;
    ijShellPairList.reserve(maxNumOfPairs);
    if (hasCalced == true) {
        return ijShellPairList;
    }

    int numOfPairs = 0;
    // I-shell Type ------------------------------------------------------
    for (int ity = state.ity; ity >= 0; --ity) {

        // J-shell Type ----------------------------------------------------
        for (int jty = state.jty; jty >= 0; --jty) {

            for (std::size_t i = std::max(0, state.i); i < shellList1[ity].size(); ++i) {
                const std::size_t index_i = shellList1[ity][i];

                for (std::size_t j = std::max(0, state.j); j < shellList2[jty].size(); ++j) {
                    const std::size_t index_j = shellList2[jty][j];

                    // cut off
                    this->countupTotal(ity, jty);
                    if (this->isCutoffUsingSchwartzInequality(orbitalInfo1, orbitalInfo2,
                                                              index_i, index_j) == true) {
                        this->countupCutoff(ity, jty);
                        continue;
                    }

                    ijShellPairList.push_back(IJShellPair(index_i,
                                                          index_j));
                    ++numOfPairs;

                    if (numOfPairs > maxNumOfPairs) {
                        // 現状の保存
                        state.ity = ity;
                        state.jty = jty;
                        state.i = i;
                        state.j = j +1; // "+1" means to start next j-shell at next time.

                        return ijShellPairList;
                    }
                }
                state.j = -1;
            }
            state.i = -1;
        }
        state.jty = MAXTYPE -1;
    }

    // finalize
    hasCalced = true;
    return ijShellPairList;
}

void DfOverlap::finalizeIntegral(TlSymmetricMatrix& rMatrix)
{
    // do nothing
    //std::cerr << "DfOverlap::finalizeIntegral(m)" << std::endl;
}

void DfOverlap::finalizeIntegral(TlVector& rVector)
{
    // do nothing
    //std::cerr << "DfOverlap::finalizeIntegral(v)" << std::endl;
}

/*****************************************************************************/
/*   Main routine of Overlap calculation for [a]                             */
/*****************************************************************************/
void DfOverlap::calcNa(const std::vector<std::size_t>& shellList, TlVector* pNa)
{
    assert(pNa != NULL);

    const double ZERO = 0.0;
    const double ONE = 1.0;
    double O[3];
    O[0] = 0.0;
    O[1] = 0.0;
    O[2] = 0.0;

    double OVP[MAX_SHELL * 2 + 1];
    std::vector<std::size_t>::const_iterator pEnd = shellList.end();
    for (std::vector<std::size_t>::const_iterator p = shellList.begin(); p != pEnd; ++p) {
        const std::size_t orbA = *p;
        const int nqA = this->pOrbitalInfo_Density_->getShellType(orbA);
        const double Ca = this->pOrbitalInfo_Density_->getCoefficient(orbA, 0);
        const double Za = this->pOrbitalInfo_Density_->getExponent(orbA, 0);
        double A[3];
        const TlPosition posA = this->pOrbitalInfo_Density_->getPosition(orbA);
        A[0] = posA.x();
        A[1] = posA.y();
        A[2] = posA.z();

        // check the index number; ic?, iz?, nu?
        for (int i = 0; i < (MAX_SHELL * 2 + 1); ++i) {
            OVP[i] = 0.0;
        }

        this->ovpqgcalc(nqA, 0, 0,
                        1, &Ca, &Za, A,
                        1, &ONE, &ZERO, A,
                        1.0, 0.0, O,
                        OVP);

        int iww = 0;
        const int iAngA = nqA * 2 +1;
        for (int iAng = 0; iAng < iAngA; ++iAng) {
            const std::size_t I = orbA + iAng;

            pNa->add(I, OVP[iww]);
            ++iww;
        }
    }
}


//   Main routine of Overlap calculation for [pq] & [gd]
void DfOverlap::spqcalc(const TlOrbitalInfoObject& orbInfo,
                        const std::vector<IJShellPair>& IJShellList, TlMatrixObject* Spq)
{
    double O[3];
    O[0] = 0.0;
    O[1] = 0.0;
    O[2] = 0.0;

    const double dCutValue = this->cutvalue;
    const int maxNumOfPairs = IJShellList.size();
#pragma omp parallel
    {
        double OVP[(MAX_SHELL * 2 + 1) *(MAX_SHELL * 2 + 1)];
#pragma omp for
        for (int pairIndex = 0; pairIndex < maxNumOfPairs; ++pairIndex) {
            const std::size_t orbA = IJShellList[pairIndex].nIShell;
            const std::size_t orbB = IJShellList[pairIndex].nJShell;
            const int nqA = orbInfo.getShellType(orbA);
            const int nqB = orbInfo.getShellType(orbB);
            
            const int npA = orbInfo.getCgtoContraction(orbA);
            double* Ca = new double[npA];
            double* Za = new double[npA];
            for (int i = 0; i < npA; ++i) {
                Ca[i] = orbInfo.getCoefficient(orbA, i);
                Za[i] = orbInfo.getExponent(orbA, i);
            }
            double A[3];
            const TlPosition posA = orbInfo.getPosition(orbA);
            A[0] = posA.x();
            A[1] = posA.y();
            A[2] = posA.z();
            
            const int npB = orbInfo.getCgtoContraction(orbB);
            double* Cb = new double[npB];
            double* Zb = new double[npB];
            for (int i = 0; i < npB; ++i) {
                Cb[i] = orbInfo.getCoefficient(orbB, i);
                Zb[i] = orbInfo.getExponent(orbB, i);
            }
            double B[3];
            const TlPosition posB = orbInfo.getPosition(orbB);
            B[0] = posB.x();
            B[1] = posB.y();
            B[2] = posB.z();
            
            for (int i = 0; i < ((MAX_SHELL * 2 + 1) *(MAX_SHELL * 2 + 1)); ++i) {
                OVP[i] = 0.0;
            }
            
            this->ovpqgcalc(nqA, nqB, 0,
                            npA, Ca, Za, A,
                            npB, Cb, Zb, B,
                            1.0, 0.0, O,
                            OVP);
            
            const int iAngA = nqA * 2 + 1;
            const int iAngB = nqB * 2 + 1;
            
#pragma omp critical (DfOverlap_spqcalc)
            {
                int iww = 0;
                for (int iAng = 0; iAng < iAngA; ++iAng) {
                    const std::size_t I = orbA + iAng;
                    
                    for (int jAng = 0; jAng < iAngB; ++jAng) {
                        const std::size_t J = orbB + jAng;
                        
                        if ((orbA != orbB) || (I >= J)) {
                            if (std::fabs(OVP[iww]) >= dCutValue) {
                                Spq->add(I, J, OVP[iww]);
                            }
                        }
                        ++iww;
                    }
                }
            }
            
            delete[] Zb;
            delete[] Cb;
            delete[] Za;
            delete[] Ca;
        }
    }
}

void DfOverlap::spqcalc(const TlOrbitalInfoObject& orbInfo1,
                        const TlOrbitalInfoObject& orbInfo2,
                        const std::vector<IJShellPair>& IJShellList, TlMatrixObject* Spq)
{
    double O[3];
    O[0] = 0.0;
    O[1] = 0.0;
    O[2] = 0.0;

    const double dCutValue = this->cutvalue;
    const int maxNumOfPairs = IJShellList.size();

#pragma omp parallel
    {
        double OVP[(MAX_SHELL * 2 + 1) *(MAX_SHELL * 2 + 1)];
#pragma omp for
        for (int pairIndex = 0; pairIndex < maxNumOfPairs; ++pairIndex) {
            std::size_t orbA = IJShellList[pairIndex].nIShell;
            std::size_t orbB = IJShellList[pairIndex].nJShell;
            
            const int nqA = orbInfo1.getShellType(orbA);
            const int npA = orbInfo1.getCgtoContraction(orbA);
            double* Ca = new double[npA];
            double* Za = new double[npA];
            for (int i = 0; i < npA; ++i) {
                Ca[i] = orbInfo1.getCoefficient(orbA, i);
                Za[i] = orbInfo1.getExponent(orbA, i);
            }
            double A[3];
            const TlPosition posA = orbInfo1.getPosition(orbA);
            A[0] = posA.x();
            A[1] = posA.y();
            A[2] = posA.z();
            
            const int nqB = orbInfo2.getShellType(orbB);
            const int npB = orbInfo2.getCgtoContraction(orbB);
            double* Cb = new double[npB];
            double* Zb = new double[npB];
            for (int i = 0; i < npB; ++i) {
                Cb[i] = orbInfo2.getCoefficient(orbB, i);
                Zb[i] = orbInfo2.getExponent(orbB, i);
            }
            double B[3];
            const TlPosition posB = orbInfo2.getPosition(orbB);
            B[0] = posB.x();
            B[1] = posB.y();
            B[2] = posB.z();
            
            for (int i = 0; i < ((MAX_SHELL * 2 + 1) *(MAX_SHELL * 2 + 1)); ++i) {
                OVP[i] = 0.0;
            }
            
            // force (nqA >= nqB)
//             bool doConvert = false;
            if (nqA >= nqB) {
                this->ovpqgcalc(nqA, nqB, 0,
                                npA, Ca, Za, A,
                                npB, Cb, Zb, B,
                                1.0, 0.0, O,
                                OVP);
            } else {
                this->ovpqgcalc(nqB, nqA, 0,
                                npB, Cb, Zb, B,
                                npA, Ca, Za, A,
                                1.0, 0.0, O,
                                OVP);

//                 doConvert = true;
            }
            
            const int iAngA = nqA * 2 + 1;
            const int iAngB = nqB * 2 + 1;
            
#pragma omp critical (DfOverlap_spqcalc_2orb)
            {
                int iww = 0;
                for (int iAng = 0; iAng < iAngA; ++iAng) {
                    const std::size_t I = orbA + iAng;
                    
                    for (int jAng = 0; jAng < iAngB; ++jAng) {
                        const std::size_t J = orbB + jAng;
                        
                        if (std::fabs(OVP[iww]) >= dCutValue) {
//                             if (doConvert == false) {
                                Spq->add(I, J, OVP[iww]);
//                             } else {
//                                 Spq->add(J, I, OVP[iww]);
//                             }
                        }
                        ++iww;
                    }
                }
            }
            
            delete[] Zb;
            delete[] Cb;
            delete[] Za;
            delete[] Ca;
        }
    }
}

void DfOverlap::ovpqgDH(const std::vector<IJShellPair>& IJShellList,
                        const TlVector& deltaMyu, TlMatrixObject* deltaH)
{
    //const std::size_t numOfAuxXC = this->numOfAuxXC_;
    const int numOfPairs = IJShellList.size();

#pragma omp parallel
    {
        double OVP[(MAX_SHELL * 2 + 1) *(MAX_SHELL * 2 + 1) *(MAX_SHELL * 2 + 1)];
#pragma omp for
        for (int pairIndex = 0; pairIndex < numOfPairs; ++pairIndex) {
            const std::size_t orbA = IJShellList[pairIndex].nIShell;
            const std::size_t orbB = IJShellList[pairIndex].nJShell;
            const int nqA = this->pOrbitalInfo_->getShellType(orbA);
            const int nqB = this->pOrbitalInfo_->getShellType(orbB);
            const int iAngA = nqA * 2 + 1;
            const int iAngB = nqB * 2 + 1;

            const int npA = this->pOrbitalInfo_->getCgtoContraction(orbA);
            double* Ca = new double[npA];
            double* Za = new double[npA];
            for (int i = 0; i < npA; ++i) {
                Ca[i] = this->pOrbitalInfo_->getCoefficient(orbA, i);
                Za[i] = this->pOrbitalInfo_->getExponent(orbA, i);
            }
            double A[3];
            const TlPosition posA = this->pOrbitalInfo_->getPosition(orbA);
            A[0] = posA.x();
            A[1] = posA.y();
            A[2] = posA.z();

            const int npB = this->pOrbitalInfo_->getCgtoContraction(orbB);
            double* Cb = new double[npB];
            double* Zb = new double[npB];
            for (int i = 0; i < npB; ++i) {
                Cb[i] = this->pOrbitalInfo_->getCoefficient(orbB, i);
                Zb[i] = this->pOrbitalInfo_->getExponent(orbB, i);
            }
            double B[3];
            const TlPosition posB = this->pOrbitalInfo_->getPosition(orbB);
            B[0] = posB.x();
            B[1] = posB.y();
            B[2] = posB.z();

            const std::size_t numOfShellType = this->shellList_Dens_.size();
            for (int shellType = numOfShellType -1; shellType >= 0; --shellType) {
                const std::size_t numOfOrbs = this->shellList_Dens_[shellType].size();
                for (std::size_t k = 0; k < numOfOrbs; ++k) {
                    const std::size_t orbC = this->shellList_Dens_[shellType][k];

                    const int nqC = this->pOrbitalInfo_XC_->getShellType(orbC);
                    const int iAngC = nqC * 2 + 1;
                    const double Cc = this->pOrbitalInfo_XC_->getCoefficient(orbC, 0);
                    const double Zc = this->pOrbitalInfo_XC_->getExponent(orbC, 0);
                    double C[3];
                    const TlPosition posC = this->pOrbitalInfo_XC_->getPosition(orbC);
                    C[0] = posC.x();
                    C[1] = posC.y();
                    C[2] = posC.z();

                    for (int i = 0; i < ((MAX_SHELL * 2 + 1) *(MAX_SHELL * 2 + 1)); ++i) {
                        OVP[i] = 0.0;
                    }

                    this->ovpqgcalc(nqA, nqB, nqC,
                                    npA, Ca, Za, A,
                                    npB, Cb, Zb, B,
                                    Cc, Zc, C,
                                    OVP);

#pragma omp critical (DfOverlap_ovpqgDH1)
                    {
                        int iww = 0;
                        for (int iAng = 0; iAng < iAngA; ++iAng) {
                            const std::size_t I = orbA + iAng;

                            for (int jAng = 0; jAng < iAngB; ++jAng) {
                                const std::size_t J = orbB + jAng;

                                double value1 = 0.0;
                                for (int kAng = 0; kAng < iAngC; ++kAng) {
                                    std::size_t K = orbC + kAng;

                                    if ((orbA != orbB) || (I >= J)) {
                                        if (std::fabs(OVP[iww]) >= cutvalue) {
                                            value1 += OVP[iww] * deltaMyu[K];
                                        }
                                    }
                                    ++iww;
                                }
                                deltaH->add(I, J, value1);
                            }
                        }
                    }

                }
            }

            delete[] Zb;
            delete[] Cb;
            delete[] Za;
            delete[] Ca;
        }
    }
}

void DfOverlap::ovpqgDH(const std::vector<IJShellPair>& IJShellList,
                        const TlVector& deltaMyu, const TlVector& deltaEps,
                        TlMatrixObject* deltaH1, TlMatrixObject* deltaH2)
{
    //const std::size_t numOfAuxXC = this->numOfAuxXC_;

    const int numOfPairs = IJShellList.size();
#pragma omp parallel
    {
        double OVP[(MAX_SHELL * 2 + 1) *(MAX_SHELL * 2 + 1) *(MAX_SHELL * 2 + 1)];
#pragma omp for
        for (int pairIndex = 0; pairIndex < numOfPairs; ++pairIndex) {
            const std::size_t orbA = IJShellList[pairIndex].nIShell;
            const std::size_t orbB = IJShellList[pairIndex].nJShell;
            const int nqA = this->pOrbitalInfo_->getShellType(orbA);
            const int nqB = this->pOrbitalInfo_->getShellType(orbB);
            const int iAngA = nqA * 2 + 1;
            const int iAngB = nqB * 2 + 1;

            const int npA = this->pOrbitalInfo_->getCgtoContraction(orbA);
            double* Ca = new double[npA];
            double* Za = new double[npA];
            for (int i = 0; i < npA; ++i) {
                Ca[i] = this->pOrbitalInfo_->getCoefficient(orbA, i);
                Za[i] = this->pOrbitalInfo_->getExponent(orbA, i);
            }
            double A[3];
            const TlPosition posA = this->pOrbitalInfo_->getPosition(orbA);
            A[0] = posA.x();
            A[1] = posA.y();
            A[2] = posA.z();

            const int npB = this->pOrbitalInfo_->getCgtoContraction(orbB);
            double* Cb = new double[npB];
            double* Zb = new double[npB];
            for (int i = 0; i < npB; ++i) {
                Cb[i] = this->pOrbitalInfo_->getCoefficient(orbB, i);
                Zb[i] = this->pOrbitalInfo_->getExponent(orbB, i);
            }
            double B[3];
            const TlPosition posB = this->pOrbitalInfo_->getPosition(orbB);
            B[0] = posB.x();
            B[1] = posB.y();
            B[2] = posB.z();


            const std::size_t numOfShellType = this->shellList_Dens_.size();
            for (int shellType = numOfShellType -1; shellType >= 0; --shellType) {
                const std::size_t numOfOrbs = this->shellList_Dens_[shellType].size();
                for (std::size_t k = 0; k < numOfOrbs; ++k) {
                    const std::size_t orbC = this->shellList_Dens_[shellType][k];

                    const int nqC = this->pOrbitalInfo_XC_->getShellType(orbC);
                    const int iAngC = nqC * 2 + 1;
                    const double Cc = this->pOrbitalInfo_XC_->getCoefficient(orbC, 0);
                    const double Zc = this->pOrbitalInfo_XC_->getExponent(orbC, 0);
                    double C[3];
                    const TlPosition posC = this->pOrbitalInfo_XC_->getPosition(orbC);
                    C[0] = posC.x();
                    C[1] = posC.y();
                    C[2] = posC.z();

                    for (int i = 0; i < ((MAX_SHELL * 2 + 1) *(MAX_SHELL * 2 + 1)); ++i) {
                        OVP[i] = 0.0;
                    }

                    this->ovpqgcalc(nqA, nqB, nqC,
                                    npA, Ca, Za, A,
                                    npB, Cb, Zb, B,
                                    Cc, Zc, C,
                                    OVP);

#pragma omp critical (DfOverlap_ovpqgDH2)
                    {
                        int iww = 0;
                        for (int iAng = 0; iAng < iAngA; ++iAng) {
                            const std::size_t I = orbA + iAng;

                            for (int jAng = 0; jAng < iAngB; ++jAng) {
                                const std::size_t J = orbB + jAng;

                                double value1 = 0.0;
                                double value2 = 0.0;
                                for (int kAng = 0; kAng < iAngC; ++kAng) {
                                    std::size_t K = orbC + kAng;

                                    if ((orbA != orbB) || (I >= J)) {
                                        if (std::fabs(OVP[iww]) >= cutvalue) {
                                            value1 += OVP[iww] * deltaMyu[K];
                                            value2 += OVP[iww] * deltaEps[K];
                                        }
                                    }
                                    ++iww;
                                }
                                deltaH1->add(I, J, value1);
                                deltaH2->add(I, J, value2);
                            }
                        }
                    }

                }
            }

            delete[] Zb;
            delete[] Cb;
            delete[] Za;
            delete[] Ca;
        }
    }
}



// Schwaltz inequality cut off routine for [pqGamma]
bool DfOverlap::isCutoffUsingSchwartzInequality(const TlOrbitalInfo& orbInfo,
                                                const std::size_t orb1, const std::size_t orb2) const
{
//     bool bAnswer = true;

//     // SS OVERLAP
//     const TlPosition A = orbInfo.getPosition(p);
//     const TlPosition B = orbInfo.getPosition(q);
//     const double absQ = A.squareDistanceFrom(B);

//     const int npA = orbInfo.getCgtoContraction(p);
//     const int npB = orbInfo.getCgtoContraction(q);

//     const double CPAI = M_PI * sqrt(M_PI);
//     double STS = 0.0;
//     for (int ipA = 0; ipA < npA; ++ipA) {
//         const double ZAW = orbInfo.getExponent(p, ipA);
//         const double coeff = CPAI * orbInfo.getCoefficient(p, ipA);

//         for (int ipB = 0; ipB < npB; ++ipB) {
//             const double ZBW = orbInfo.getExponent(q, ipB);
//             const double ZP = ZAW + ZBW;
//             const double ZPIW = 1.0 / ZP;
//             const double coefB = orbInfo.getCoefficient(q, ipB);
//             STS += coeff * ZPIW * sqrt(ZPIW) * exp(-absQ * ZAW * ZBW * ZPIW) * coefB;
//         }
//     }

//     if (std::fabs(STS) >= this->cutvalue) {
//         bAnswer = false;
//     }

//     return bAnswer;
    return this->isCutoffUsingSchwartzInequality(orbInfo, orbInfo, orb1, orb2);
}


bool DfOverlap::isCutoffUsingSchwartzInequality(const TlOrbitalInfo& orbInfo1,
                                                const TlOrbitalInfo& orbInfo2,
                                                const std::size_t orb1,
                                                const std::size_t orb2) const
{
    bool bAnswer = true;

    // SS OVERLAP
    const TlPosition A = orbInfo1.getPosition(orb1);
    const TlPosition B = orbInfo2.getPosition(orb2);
    const double absQ = A.squareDistanceFrom(B);

    const int npA = orbInfo1.getCgtoContraction(orb1);
    const int npB = orbInfo2.getCgtoContraction(orb2);

    const double CPAI = M_PI * sqrt(M_PI);
    double STS = 0.0;
    for (int ipA = 0; ipA < npA; ++ipA) {
        const double ZAW = orbInfo1.getExponent(orb1, ipA);
        const double coeff = CPAI * orbInfo1.getCoefficient(orb1, ipA);

        for (int ipB = 0; ipB < npB; ++ipB) {
            const double ZBW = orbInfo2.getExponent(orb2, ipB);
            const double ZP = ZAW + ZBW;
            const double ZPIW = 1.0 / ZP;
            const double coefB = orbInfo2.getCoefficient(orb2, ipB);
            STS += coeff * ZPIW * sqrt(ZPIW) * exp(-absQ * ZAW * ZBW * ZPIW) * coefB;
        }
    }

    if (std::fabs(STS) >= this->cutvalue) {
        bAnswer = false;
    }

    return bAnswer;
}


/*****************************************************************************/
/*   Overlap calculation for [pqGamma] ( [ABC] )                 */
/*  Input:  nqA ; Quntum number of A orbital. 0, 1, 2, ...       */
/*      nqB ; Quntum number of B orbital. 0, 1, 2, ...       */
/*      nqC ; Quntum number of C orbital. 0, 1, 2, ...       */
/*      npA ; Number of primitives in A              */
/*      npB ; Number of primitives in B              */
/*      CA  ; Contraction coefficients for A             */
/*      CB  ; Contraction coefficients for B             */
/*      CC  ; Contraction coefficient for C              */
/*      ZA  ; Orbital exponents for A                */
/*      ZB  ; Orbital exponents for B                */
/*      ZC  ; Orbital exponent for C                 */
/*      A,B,C   ; Nuclear coordinates                    */
/*****************************************************************************/
void DfOverlap::ovpqgcalc(const int nqA, const int nqB, const int nqC,
                          const int npA, const double* CA, const double* ZA, const double* A,
                          const int npB, const double* CB, const double* ZB, const double* B,
                          const double CC, const double ZC, const double* C,
                          double* OVP)
{
    assert((nqA < 3) && (nqB < 3) && (nqC < 3));
    assert(nqA >= nqB);
    assert((npA <= 10) && (npB <= 10));

    const int iAngA = nqA * 2 + 1;
    const int iAngB = nqB * 2 + 1;
    const int iAngC = nqC * 2 + 1;
    const int nAng  = iAngA*iAngB*iAngC;

    double absQ = 0.0;
    for (int i = 0; i < 3; ++i) {
        const double ABW   = A[i] - B[i];
        absQ += (ABW * ABW);
    }

    // allocate ERP(work memory)
    const int ERP_R_SIZE = 101;
    const int ERP_PQ_SIZE = 126;
    double** ERP = new double*[ERP_R_SIZE];
    for (int i = 0; i < ERP_R_SIZE; ++i) {
        ERP[i] = new double[ERP_PQ_SIZE];
    }

    const int numOfPQContracts = npA * npB;
    //double P[3][101], PA[3][101], PB[3][101];
    double** P = new double*[3];
    double** PA = new double*[3];
    double** PB = new double*[3];
    for (int i = 0; i < 3; ++i) {
        P[i] = new double[numOfPQContracts];
        PA[i] = new double[numOfPQContracts];
        PB[i] = new double[numOfPQContracts];
    }
    double* Zp = new double[numOfPQContracts];
    double* HP = new double[numOfPQContracts];

    int ip = 0;
    for (int ipA = 0; ipA < npA; ++ipA) {
        for (int ipB = 0; ipB < npB; ++ipB) {
            const double ZAW = ZA[ipA];
            const double ZBW = ZB[ipB];
            Zp[ip]  = ZAW + ZBW;
            const double ZPIW = 1.0 / Zp[ip];
            HP[ip]  = exp(-absQ * ZAW * ZBW * ZPIW) * CA[ipA] * CB[ipB];
            for (int i = 0; i < 3; ++i) {
                P [i][ip] = (ZAW*A[i] + ZBW*B[i])*ZPIW;
                PA[i][ip] = A[i];
                PB[i][ip] = B[i];
            }
            ++ip;
        }
    }

    const int np = npA * npB;
    const int judge = nqA * 16 + nqB * 4 + nqC;
    switch (judge) {
    case 0:
        this->ovpSSS(np, P[0], P[1], P[2], ZC, C[0], C[1], C[2], CC, Zp, HP, ERP);
        break;

    case 1:
        this->ovpSSP(np, P[0], P[1], P[2], ZC, C[0], C[1], C[2], CC, Zp, HP, ERP);
        break;

    case 2:
        this->ovpSSD(np, P[0], P[1], P[2], ZC, C[0], C[1], C[2], CC, Zp, HP, ERP);
        break;

    case 16:
        this->ovpPSS(np, P[0], P[1], P[2], PA[0], PA[1], PA[2], ZC, C[0], C[1], C[2], CC, Zp, HP, ERP);
        break;

    case 17:
        this->ovpPSP(np, P[0], P[1], P[2], PA[0], PA[1], PA[2], ZC, C[0], C[1], C[2], CC, Zp, HP, ERP);
        break;

    case 18:
        this->ovpPSD(np, P[0], P[1], P[2], PA[0], PA[1], PA[2], ZC, C[0], C[1], C[2], CC, Zp, HP, ERP);
        break;

    case 20:
        this->ovpPPS(np, P[0], P[1], P[2], PA[0], PA[1], PA[2], PB[0], PB[1], PB[2], ZC, C[0], C[1], C[2], CC, Zp, HP, ERP);
        break;

    case 21:
        this->ovpPPP(np, P[0], P[1], P[2], PA[0], PA[1], PA[2], PB[0], PB[1], PB[2], ZC, C[0], C[1], C[2], CC, Zp, HP, ERP);
        break;

    case 22:
        this->ovpPPD(np, P[0], P[1], P[2], PA[0], PA[1], PA[2], PB[0], PB[1], PB[2], ZC, C[0], C[1], C[2], CC, Zp, HP, ERP);
        break;

    case 32:
        this->ovpDSS(np, P[0], P[1], P[2], PA[0], PA[1], PA[2], ZC, C[0], C[1], C[2], CC, Zp, HP, ERP);
        break;

    case 33:
        this->ovpDSP(np, P[0], P[1], P[2], PA[0], PA[1], PA[2], ZC, C[0], C[1], C[2], CC, Zp, HP, ERP);
        break;

    case 34:
        this->ovpDSD(np, P[0], P[1], P[2], PA[0], PA[1], PA[2], ZC, C[0], C[1], C[2], CC, Zp, HP, ERP);
        break;

    case 36:
        this->ovpDPS(np, P[0], P[1], P[2], PA[0], PA[1], PA[2], PB[0], PB[1], PB[2], ZC, C[0], C[1], C[2], CC, Zp, HP, ERP);
        break;

    case 37:
        this->ovpDPP(np, P[0], P[1], P[2], PA[0], PA[1], PA[2], PB[0], PB[1], PB[2], ZC, C[0], C[1], C[2], CC, Zp, HP, ERP);
        break;

    case 38:
        this->ovpDPD(np, P[0], P[1], P[2], PA[0], PA[1], PA[2], PB[0], PB[1], PB[2], ZC, C[0], C[1], C[2], CC, Zp, HP, ERP);
        break;

    case 40:
        this->ovpDDS(np, P[0], P[1], P[2], PA[0], PA[1], PA[2], PB[0], PB[1], PB[2], ZC, C[0], C[1], C[2], CC, Zp, HP, ERP);
        break;

    case 41:
        this->ovpDDP(np, P[0], P[1], P[2], PA[0], PA[1], PA[2], PB[0], PB[1], PB[2], ZC, C[0], C[1], C[2], CC, Zp, HP, ERP);
        break;

    case 42:
        this->ovpDDD(np, P[0], P[1], P[2], PA[0], PA[1], PA[2], PB[0], PB[1], PB[2], ZC, C[0], C[1], C[2], CC, Zp, HP, ERP);
        break;

    default:
        abort();
        break;
    }

    // Contraction loop
    for (int i = 0; i < nAng; ++i) {
        OVP[i] = 0.0;
    }

    for (int i = 0; i < nAng; ++i) {
        for (int ip = 0; ip < np; ++ip) {
            OVP[i] += ERP[ip][i];
        }
    }

    // delete
    delete[] Zp;
    delete[] HP;

    for (int i = 0; i < 3; ++i) {
        delete[] P[i];
        delete[] PA[i];
        delete[] PB[i];
    }
    delete[] P;
    delete[] PA;
    delete[] PB;

    for (int i = 0; i < ERP_R_SIZE; ++i) {
        delete[] ERP[i];
    }
    delete[] ERP;
}


/*****************************************************************************/
/*   make table for OVP [pqGamma]                        */
/*                                       */
/*  natom  = total number of atoms                       */
/*  ncgto  =         orbital CGTOs                   */
/*  nexgto =         exchange function GTOs              */
/*  norbf  =         orbital basis functions             */
/*  nauxxc =         exchange functions              */
/*****************************************************************************/
void DfOverlap::makeTable()
{
//     const std::size_t numOfAOs = this->m_nNumOfAOs;
//     this->shellList_.clear();
//     this->shellList_.resize(3); // s, p, d support
//     {
//         std::size_t shell = 0;
//         while (shell < numOfAOs) {
//             const int shellType = this->pOrbitalInfo_->getShellType(shell);
//             this->shellList_[shellType].push_back(shell);
//             shell += shellType * 2 + 1;
//         }
//     }
    this->shellList_ = this->makeShellList(*(this->pOrbitalInfo_));

//     const std::size_t numOfAuxDens = this->m_nNumOfAux;
//     this->shellList_Dens_.clear();
//     this->shellList_Dens_.resize(3); // s, p, d support
//     {
//         std::size_t shell = 0;
//         while (shell < numOfAuxDens) {
//             const int shellType = this->pOrbitalInfo_Density_->getShellType(shell);
//             this->shellList_Dens_[shellType].push_back(shell);
//             shell += shellType * 2 + 1;
//         }
//     }
    this->shellList_Dens_ = this->makeShellList(*(this->pOrbitalInfo_Density_));

//     const std::size_t numOfAuxXC = this->numOfAuxXC_;
//     this->shellList_XC_.clear();
//     this->shellList_XC_.resize(3); // s, p, d support
//     {
//         std::size_t shell = 0;
//         while (shell < numOfAuxXC) {
//             const int shellType = this->pOrbitalInfo_XC_->getShellType(shell);
//             this->shellList_XC_[shellType].push_back(shell);
//             shell += shellType * 2 + 1;
//         }
//     }
    this->shellList_XC_ = this->makeShellList(*(this->pOrbitalInfo_XC_));
}

DfOverlap::ShellListType DfOverlap::makeShellList(const TlOrbitalInfoObject& orbInfo)
{
    const std::size_t numOfOrbitals = orbInfo.getNumOfOrbitals();
    const int MaxShellType = 3; // s, p, d

    ShellListType shellList(MaxShellType);

    std::size_t index = 0;
    while (index < numOfOrbitals) {
        const int shellType = orbInfo.getShellType(index);
        shellList[shellType].push_back(index);
        index += shellType * 2 + 1;
    }
    
    return shellList;
}

/*****************************************************************************/
/*  Primitive 3-center OVP evaluation [ ABC ] for type [ SSS ]           */
/*  np      ; number of first primitive pair data            */
/*  Zp      ; (Alpha+Beta) from orbital exponents Alpha, Beta    */
/*  px, py, pz  ; (Alpha*ax * Beta*bx) / (Alpha+Beta) , ...          */
/*  HP      ; Zeta**(-3/2) * exp(-(Alpha*Beta/Zeta) * (A-B)**2   */
/*  Gamma       ; orbital exponent Gamma for last GTO            */
/*  cx, cy, cz  ; coordinates of last GTO C              */
/*  cc      ; coefficient of last GTO C              */
/*  ERP     ; output primitive OVP                   */
/*****************************************************************************/
void DfOverlap::ovpSSS(const int np,
                       const double* px, const double* py, const double* pz, const double Gamma,
                       const double cx, const double cy, const double cz, const double cc,
                       const double* Zp, const double* HP, double** ERP)
{
    const double hc = M_PI*sqrt(M_PI)*cc;
    for (int i = 0; i < np; ++i) {
        const double RZG = 1.0 / (Zp[i] + Gamma);
        const double Roo = Zp[i]*Gamma*RZG;
        const double pcx = px[i] - cx;
        const double pcy = py[i] - cy;
        const double pcz = pz[i] - cz;
        const double t = Roo * (pcx*pcx + pcy*pcy + pcz*pcz);

        ERP[i][0] = hc * HP[i] * RZG * sqrt(RZG) * exp(-t);
    }
}

/*****************************************************************************/
/*  Primitive 3-center OVP evaluation [ ABC ] for type [ PSS ]           */
/*  np      ; number of first primitive pair data            */
/*  Zp      ; (Alpha+Beta) from orbital exponents Alpha, Beta    */
/*  px, py, pz  ; (Alpha*ax * Beta*bx) / (Alpha+Beta) , ...          */
/*  pax, pay, paz   ; ax, ay, az                         */
/*  HP      ; Zeta**(-3/2) * exp(-(Alpha*Beta/Zeta) * (A-B)**2   */
/*  Gamma       ; orbital exponent Gamma for last GTO            */
/*  cx, cy, cz  ; coordinates of last GTO C              */
/*  cc      ; coefficient of last GTO C              */
/*  ERP     ; output primitive OVP                   */
/*              ERP[0,1,2,...,np-1][0] for [ PxSS ]      */
/*              ERP[0,1,2,...,np-1][1] for [ PySS ]      */
/*              ERP[0,1,2,...,np-1][2] for [ PzSS ]      */
/*****************************************************************************/
void DfOverlap::ovpPSS(const int np,
                       const double* px, const double* py, const double* pz,
                       const double* pax, const double* pay, const double* paz,
                       const double Gamma,
                       const double cx, const double cy, const double cz,
                       const double cc,
                       const double* Zp, const double* HP, double** ERP)
{
    const double hc = M_PI*sqrt(M_PI)*cc;
    for (int i = 0; i < np; ++i) {
        const double RZG = 1.0 / (Zp[i] + Gamma);
        const double Roo = Zp[i]*Gamma*RZG;
        const double pcx = px[i] - cx;
        const double pcy = py[i] - cy;
        const double pcz = pz[i] - cz;
        const double t = Roo * (pcx*pcx + pcy*pcy + pcz*pcz);
        const double sss = hc * HP[i] * RZG * sqrt(RZG) * exp(-t);

        const double Gx = (Zp[i]*px[i] + Gamma*cx) * RZG;
        const double Gy = (Zp[i]*py[i] + Gamma*cy) * RZG;
        const double Gz = (Zp[i]*pz[i] + Gamma*cz) * RZG;
        const double Gax = Gx - pax[i];
        const double Gay = Gy - pay[i];
        const double Gaz = Gz - paz[i];

        ERP[i][0] = Gax*sss;
        ERP[i][1] = Gay*sss;
        ERP[i][2] = Gaz*sss;
    }
}

/*****************************************************************************/
/*                                       */
/*  Primitive 3-center OVP evaluation [ ABC ] for type [ SSP ]           */
/*                                       */
/*  np      ; number of first primitive pair data            */
/*  Zp      ; (Alpha+Beta) from orbital exponents Alpha, Beta    */
/*  px, py, pz  ; (Alpha*ax * Beta*bx) / (Alpha+Beta) , ...          */
/*  HP      ; Zeta**(-3/2) * exp(-(Alpha*Beta/Zeta) * (A-B)**2   */
/*  Gamma       ; orbital exponent Gamma for last GTO            */
/*  cx, cy, cz  ; coordinates of last GTO C              */
/*  cc      ; coefficient of last GTO C              */
/*  ERP     ; output primitive OVP                   */
/*              ERP[0,1,2,...,np-1][0] for [ SSPx ]      */
/*              ERP[0,1,2,...,np-1][1] for [ SSPy ]      */
/*              ERP[0,1,2,...,np-1][2] for [ SSPz ]      */
/*                                       */
/*****************************************************************************/
void DfOverlap::ovpSSP(const int np,
                       const double* px, const double* py, const double* pz,
                       const double Gamma,
                       const double cx, const double cy, const double cz,
                       const double cc,
                       const double* Zp, const double* HP, double** ERP)
{
    const double hc = M_PI*sqrt(M_PI)*cc;
    for (int i = 0; i < np; ++i) {
        const double RZG = 1.0 / (Zp[i] + Gamma);
        const double Roo = Zp[i]*Gamma*RZG;
        const double pcx = px[i] - cx;
        const double pcy = py[i] - cy;
        const double pcz = pz[i] - cz;
        const double t   = Roo * (pcx*pcx + pcy*pcy + pcz*pcz);
        const double sss = hc * HP[i] * RZG * sqrt(RZG) * exp(-t);

        const double Gx  = (Zp[i]*px[i] + Gamma*cx) * RZG;
        const double Gy  = (Zp[i]*py[i] + Gamma*cy) * RZG;
        const double Gz  = (Zp[i]*pz[i] + Gamma*cz) * RZG;
        const double Gcx = Gx - cx;
        const double Gcy = Gy - cy;
        const double Gcz = Gz - cz;

        ERP[i][0] = Gcx*sss;
        ERP[i][1] = Gcy*sss;
        ERP[i][2] = Gcz*sss;
    }
}

/*****************************************************************************/
/*  Primitive 3-center OVP evaluation [ ABC ] for type [ PPS ]           */
/*  np      ; number of first primitive pair data            */
/*  Zp      ; (Alpha+Beta) from orbital exponents Alpha, Beta    */
/*  px, py, pz  ; (Alpha*ax * Beta*bx) / (Alpha+Beta) , ...          */
/*  pax, pay, paz   ; ax, ay, az                         */
/*  pbx, pby, pbz   ; bx, by, bz                         */
/*  HP      ; Zeta**(-3/2) * exp(-(Alpha*Beta/Zeta) * (A-B)**2   */
/*  Gamma       ; orbital exponent Gamma for last GTO            */
/*  cx, cy, cz  ; coordinates of last GTO C              */
/*  cc      ; coefficient of last GTO C              */
/*  ERP     ; output primitive OVP                   */
/*              ERP[0,1,2,...,np-1][0] for [ PxPxS ]         */
/*              ERP[0,1,2,...,np-1][1] for [ PxPyS ]         */
/*              ERP[0,1,2,...,np-1][2] for [ PxPzS ]         */
/*              ERP[0,1,2,...,np-1][3] for [ PyPxS ]         */
/*              ERP[0,1,2,...,np-1][4] for [ PyPyS ]         */
/*              ERP[0,1,2,...,np-1][5] for [ PyPzS ]         */
/*              ERP[0,1,2,...,np-1][6] for [ PzPxS ]         */
/*              ERP[0,1,2,...,np-1][7] for [ PzPyS ]         */
/*              ERP[0,1,2,...,np-1][8] for [ PzPzS ]         */
/*****************************************************************************/
void DfOverlap::ovpPPS(const int np,
                       const double* px, const double* py, const double* pz,
                       const double* pax, const double* pay, const double* paz,
                       const double* pbx, const double* pby, const double* pbz,
                       const double Gamma,
                       const double cx, const double cy, const double cz,
                       const double cc,
                       const double* Zp, const double* HP, double** ERP)
{
    const double hc  = M_PI*sqrt(M_PI)*cc;
    for (int i = 0; i < np; ++i) {
        const double RZG = 1.0 / (Zp[i] + Gamma);
        const double Roo = Zp[i]*Gamma*RZG;
        const double pcx = px[i] - cx;
        const double pcy = py[i] - cy;
        const double pcz = pz[i] - cz;
        const double t = Roo * (pcx*pcx + pcy*pcy + pcz*pcz);
        const double sss = hc * HP[i] * RZG * sqrt(RZG) * exp(-t);

        const double Gx  = (Zp[i]*px[i] + Gamma*cx) * RZG;
        const double Gy  = (Zp[i]*py[i] + Gamma*cy) * RZG;
        const double Gz  = (Zp[i]*pz[i] + Gamma*cz) * RZG;
        const double Gax = Gx - pax[i];
        const double Gay = Gy - pay[i];
        const double Gaz = Gz - paz[i];
        const double Gbx = Gx - pbx[i];
        const double Gby = Gy - pby[i];
        const double Gbz = Gz - pbz[i];
        // [ PSS ]
        const double RZG5 = 0.5*RZG;
        const double xss  = Gax*sss;
        const double yss  = Gay*sss;
        const double zss  = Gaz*sss;
        // [ PPS ]
        const double RZGS = RZG5*sss;
        ERP[i][0] = Gbx*xss + RZGS;
        ERP[i][1] = Gby*xss;
        ERP[i][2] = Gbz*xss;
        ERP[i][3] = Gbx*yss;
        ERP[i][4] = Gby*yss + RZGS;
        ERP[i][5] = Gbz*yss;
        ERP[i][6] = Gbx*zss;
        ERP[i][7] = Gby*zss;
        ERP[i][8] = Gbz*zss + RZGS;
    }
}

/*****************************************************************************/
/*  Primitive 3-center OVP evaluation [ ABC ] for type [ PSP ]           */
/*  np      ; number of first primitive pair data            */
/*  Zp      ; (Alpha+Beta) from orbital exponents Alpha, Beta    */
/*  px, py, pz  ; (Alpha*ax * Beta*bx) / (Alpha+Beta) , ...          */
/*  pax, pay, paz   ; ax, ay, az                         */
/*  HP      ; Zeta**(-3/2) * exp(-(Alpha*Beta/Zeta) * (A-B)**2   */
/*  Gamma       ; orbital exponent Gamma for last GTO            */
/*  cx, cy, cz  ; coordinates of last GTO C              */
/*  cc      ; coefficient of last GTO C              */
/*  ERP     ; output primitive OVP                   */
/*              ERP[0,1,2,...,np-1][0] for [ PxSPx ]         */
/*              ERP[0,1,2,...,np-1][1] for [ PySPx ]         */
/*              ERP[0,1,2,...,np-1][2] for [ PzSPx ]         */
/*              ERP[0,1,2,...,np-1][3] for [ PxSPy ]         */
/*              ERP[0,1,2,...,np-1][4] for [ PySPy ]         */
/*              ERP[0,1,2,...,np-1][5] for [ PzSPy ]         */
/*              ERP[0,1,2,...,np-1][6] for [ PxSPz ]         */
/*              ERP[0,1,2,...,np-1][7] for [ PySPz ]         */
/*              ERP[0,1,2,...,np-1][8] for [ PzSPz ]         */
/*****************************************************************************/
void DfOverlap::ovpPSP(const int np,
                       const double* px, const double* py, const double* pz,
                       const double* pax, const double* pay, const double* paz,
                       const double Gamma,
                       const double cx, const double cy, const double cz,
                       const double cc,
                       const double* Zp, const double* HP, double** ERP)
{
    const double hc = M_PI*sqrt(M_PI)*cc;
    for (int i = 0; i < np; ++i) {
        const double RZG = 1.0 / (Zp[i] + Gamma);
        const double Roo = Zp[i]*Gamma*RZG;
        const double pcx = px[i] - cx;
        const double pcy = py[i] - cy;
        const double pcz = pz[i] - cz;
        const double t   = Roo * (pcx*pcx + pcy*pcy + pcz*pcz);
        const double sss = hc * HP[i] * RZG * sqrt(RZG) * exp(-t);

        const double Gx  = (Zp[i]*px[i] + Gamma*cx) * RZG;
        const double Gy  = (Zp[i]*py[i] + Gamma*cy) * RZG;
        const double Gz  = (Zp[i]*pz[i] + Gamma*cz) * RZG;
        const double Gax = Gx - pax[i];
        const double Gay = Gy - pay[i];
        const double Gaz = Gz - paz[i];
        const double Gcx = Gx - cx;
        const double Gcy = Gy - cy;
        const double Gcz = Gz - cz;
        // [ PSS ]
        const double RZG5 = 0.5*RZG;
        const double xss  = Gax*sss;
        const double yss  = Gay*sss;
        const double zss  = Gaz*sss;
        // [ PSP ]
        const double RZGS = RZG5*sss;
        ERP[i][0] = Gcx*xss + RZGS;
        ERP[i][1] = Gcy*xss;
        ERP[i][2] = Gcz*xss;
        ERP[i][3] = Gcx*yss;
        ERP[i][4] = Gcy*yss + RZGS;
        ERP[i][5] = Gcz*yss;
        ERP[i][6] = Gcx*zss;
        ERP[i][7] = Gcy*zss;
        ERP[i][8] = Gcz*zss + RZGS;
    }
}

/*****************************************************************************/
/*  Primitive 3-center OVP evaluation [ ABC ] for type [ PPP ]           */
/*  np      ; number of first primitive pair data            */
/*  Zp      ; (Alpha+Beta) from orbital exponents Alpha, Beta    */
/*  px, py, pz  ; (Alpha*ax * Beta*bx) / (Alpha+Beta) , ...          */
/*  pax, pay, paz   ; ax, ay, az                         */
/*  pbx, pby, pbz   ; bx, by, bz                         */
/*  HP      ; Zeta**(-3/2) * exp(-(Alpha*Beta/Zeta) * (A-B)**2   */
/*  Gamma       ; orbital exponent Gamma for last GTO            */
/*  cx, cy, cz  ; coordinates of last GTO C              */
/*  cc      ; coefficient of last GTO C              */
/*  ERP     ; output primitive OVP                   */
/*              ERP[0,1,2,...,np-1][0] for [ PxPxPx ]        */
/*              ERP[0,1,2,...,np-1][1] for [ PxPyPx ]        */
/*              ERP[0,1,2,...,np-1][2] for [ PxPzPx ]        */
/*              ERP[0,1,2,...,np-1][3] for [ PyPxPx ]        */
/*              ERP[0,1,2,...,np-1][4] for [ PyPyPx ]        */
/*              ERP[0,1,2,...,np-1][5] for [ PyPzPx ]        */
/*              ERP[0,1,2,...,np-1][6] for [ PzPxPx ]        */
/*              ERP[0,1,2,...,np-1][7] for [ PzPyPx ]        */
/*              ERP[0,1,2,...,np-1][8] for [ PzPzPx ]        */
/*                  .............                */
/*****************************************************************************/
void DfOverlap::ovpPPP(const int np,
                       const double* px, const double* py, const double* pz,
                       const double* pax, const double* pay, const double* paz,
                       const double* pbx, const double* pby, const double* pbz,
                       const double Gamma,
                       const double cx, const double cy, const double cz,
                       const double cc,
                       const double* Zp, const double* HP, double** ERP)
{
    const double hc = M_PI*sqrt(M_PI)*cc;
    for (int i = 0; i < np; ++i) {
        const double RZG = 1.0 / (Zp[i] + Gamma);
        const double Roo = Zp[i]*Gamma*RZG;
        const double pcx = px[i] - cx;
        const double pcy = py[i] - cy;
        const double pcz = pz[i] - cz;
        const double t   = Roo * (pcx*pcx + pcy*pcy + pcz*pcz);
        const double sss = hc * HP[i] * RZG * sqrt(RZG) * exp(-t);

        const double Gx  = (Zp[i]*px[i] + Gamma*cx) * RZG;
        const double Gy  = (Zp[i]*py[i] + Gamma*cy) * RZG;
        const double Gz  = (Zp[i]*pz[i] + Gamma*cz) * RZG;
        const double Gax = Gx - pax[i];
        const double Gay = Gy - pay[i];
        const double Gaz = Gz - paz[i];
        const double Gbx = Gx - pbx[i];
        const double Gby = Gy - pby[i];
        const double Gbz = Gz - pbz[i];
        const double Gcx = Gx - cx;
        const double Gcy = Gy - cy;
        const double Gcz = Gz - cz;
        // [ PSS ]
        const double RZG5 = 0.5*RZG;
        const double xss  = Gax*sss;
        const double yss  = Gay*sss;
        const double zss  = Gaz*sss;
        // [ PSP ]
        const double sxs  = Gbx*sss;
        const double sys  = Gby*sss;
        const double szs  = Gbz*sss;
        // [ PPS ]
        const double RZGS = RZG5*sss;
        const double xxs  = Gbx*xss + RZGS;
        const double xys  = Gby*xss;
        const double xzs  = Gbz*xss;
        const double yxs  = Gbx*yss;
        const double yys  = Gby*yss + RZGS;
        const double yzs  = Gbz*yss;
        const double zxs  = Gbx*zss;
        const double zys  = Gby*zss;
        const double zzs  = Gbz*zss + RZGS;
        // [ PPP ]
        const double Rxss = RZG5*xss;
        const double Ryss = RZG5*yss;
        const double Rzss = RZG5*zss;
        const double Rsxs = RZG5*sxs;
        const double Rsys = RZG5*sys;
        const double Rszs = RZG5*szs;

        ERP[i][ 0] = Gcx*xxs + Rsxs + Rxss;
        ERP[i][ 1] = Gcy*xxs;
        ERP[i][ 2] = Gcz*xxs;
        ERP[i][ 3] = Gcx*xys + Rsys;
        ERP[i][ 4] = Gcy*xys        + Rxss;
        ERP[i][ 5] = Gcz*xys;
        ERP[i][ 6] = Gcx*xzs + Rszs;
        ERP[i][ 7] = Gcy*xzs;
        ERP[i][ 8] = Gcz*xzs        + Rxss;
        ERP[i][ 9] = Gcx*yxs        + Ryss;
        ERP[i][10] = Gcy*yxs + Rsxs;
        ERP[i][11] = Gcz*yxs;
        ERP[i][12] = Gcx*yys;
        ERP[i][13] = Gcy*yys + Rsys + Ryss;
        ERP[i][14] = Gcz*yys;
        ERP[i][15] = Gcx*yzs;
        ERP[i][16] = Gcy*yzs + Rszs;
        ERP[i][17] = Gcz*yzs        + Ryss;
        ERP[i][18] = Gcx*zxs        + Rzss;
        ERP[i][19] = Gcy*zxs;
        ERP[i][20] = Gcz*zxs + Rsxs;
        ERP[i][21] = Gcx*zys;
        ERP[i][22] = Gcy*zys        + Rzss;
        ERP[i][23] = Gcz*zys + Rsys;
        ERP[i][24] = Gcx*zzs;
        ERP[i][25] = Gcy*zzs;
        ERP[i][26] = Gcz*zzs + Rszs + Rzss;
    }
}

/*****************************************************************************/
/*  Primitive 3-center OVP evaluation [ ABC ] for type [ DSS ]           */
/*  np      ; number of first primitive pair data            */
/*  Zp      ; (Alpha+Beta) from orbital exponents Alpha, Beta    */
/*  px, py, pz  ; (Alpha*ax * Beta*bx) / (Alpha+Beta) , ...          */
/*  pax, pay, paz   ; ax, ay, az                         */
/*  HP      ; Zeta**(-3/2) * exp(-(Alpha*Beta/Zeta) * (A-B)**2   */
/*  Gamma       ; orbital exponent Gamma for last GTO            */
/*  cx, cy, cz  ; coordinates of last GTO C              */
/*  cc      ; coefficient of last GTO C              */
/*  ERP     ; output primitive OVP                   */
/*              ERP[0,1,2,...,np-1][0] for [ DxySS ]         */
/*              ERP[0,1,2,...,np-1][1] for [ DxzSS ]         */
/*              ERP[0,1,2,...,np-1][2] for [ DyzSS ]         */
/*              ERP[0,1,2,...,np-1][3] for [ Dxx-yySS ]      */
/*              ERP[0,1,2,...,np-1][4] for [ D3zz-rrSS ]     */
/*****************************************************************************/
void DfOverlap::ovpDSS(const int np,
                       const double* px, const double* py, const double* pz,
                       const double* pax, const double* pay, const double* paz,
                       const double Gamma,
                       const double cx, const double cy, const double cz,
                       const double cc,
                       const double* Zp, const double* HP, double** ERP)
{
    const double hc  = M_PI*sqrt(M_PI)*cc;
    for (int i = 0; i < np; ++i) {
        const double RZG = 1.0 / (Zp[i] + Gamma);
        const double Roo = Zp[i]*Gamma*RZG;
        const double pcx = px[i] - cx;
        const double pcy = py[i] - cy;
        const double pcz = pz[i] - cz;
        const double t   = Roo * (pcx*pcx + pcy*pcy + pcz*pcz);
        const double sss = hc * HP[i] * RZG * sqrt(RZG) * exp(-t);

        const double Gx  = (Zp[i]*px[i] + Gamma*cx) * RZG;
        const double Gy  = (Zp[i]*py[i] + Gamma*cy) * RZG;
        const double Gz  = (Zp[i]*pz[i] + Gamma*cz) * RZG;
        const double Gax = Gx - pax[i];
        const double Gay = Gy - pay[i];
        const double Gaz = Gz - paz[i];
        // [ PSS ]
        const double xss = Gax*sss;
        const double yss = Gay*sss;
        const double zss = Gaz*sss;
        // [ DSS ]  xy, xz, yz, xx-yy, 3zz-rr
        const double xxssww = Gax*xss;
        const double yyssww = Gay*yss;

        ERP[i][0] = Gax*yss;
        ERP[i][1] = Gax*zss;
        ERP[i][2] = Gay*zss;
        ERP[i][3] = 0.5*(xxssww - yyssww);
        ERP[i][4] = SQR3I*(Gaz*zss
                           - 0.5*(xxssww + yyssww));
    }
}

/*****************************************************************************/
/*  Primitive 3-center OVP evaluation [ ABC ] for type [ SSD ]           */
/*  np      ; number of first primitive pair data            */
/*  Zp      ; (Alpha+Beta) from orbital exponents Alpha, Beta    */
/*  px, py, pz  ; (Alpha*ax * Beta*bx) / (Alpha+Beta) , ...          */
/*  HP      ; Zeta**(-3/2) * exp(-(Alpha*Beta/Zeta) * (A-B)**2   */
/*  Gamma       ; orbital exponent Gamma for last GTO            */
/*  cx, cy, cz  ; coordinates of last GTO C              */
/*  cc      ; coefficient of last GTO C              */
/*  ERP     ; output primitive OVP                   */
/*              ERP[0,1,2,...,np-1][0] for [ SSDxy ]         */
/*              ERP[0,1,2,...,np-1][1] for [ SSDxz ]         */
/*              ERP[0,1,2,...,np-1][2] for [ SSDyz ]         */
/*              ERP[0,1,2,...,np-1][3] for [ SSDxx-yy ]      */
/*              ERP[0,1,2,...,np-1][4] for [ SSD3zz-rr ]     */
/*****************************************************************************/
void DfOverlap::ovpSSD(const int np,
                       const double* px, const double* py, const double* pz,
                       const double Gamma,
                       const double cx, const double cy, const double cz,
                       const double cc,
                       const double* Zp, const double* HP, double** ERP)
{
    const double hc  = M_PI*sqrt(M_PI)*cc;
    for (int i = 0; i < np; ++i) {
        const double RZG = 1.0 / (Zp[i] + Gamma);
        const double Roo = Zp[i]*Gamma*RZG;
        const double pcx = px[i] - cx;
        const double pcy = py[i] - cy;
        const double pcz = pz[i] - cz;
        const double t   = Roo * (pcx*pcx + pcy*pcy + pcz*pcz);
        const double sss = hc * HP[i] * RZG * sqrt(RZG) * exp(-t);

        const double Gx  = (Zp[i]*px[i] + Gamma*cx) * RZG;
        const double Gy  = (Zp[i]*py[i] + Gamma*cy) * RZG;
        const double Gz  = (Zp[i]*pz[i] + Gamma*cz) * RZG;
        const double Gcx = Gx - cx;
        const double Gcy = Gy - cy;
        const double Gcz = Gz - cz;
        // [ SSP ]
        const double ssx = Gcx*sss;
        const double ssy = Gcy*sss;
        const double ssz = Gcz*sss;
        // [ SSD ]  xy, xz, yz, xx-yy, 3zz-rr
        const double xssxww = Gcx*ssx;
        const double yssyww = Gcy*ssy;

        ERP[i][0] = Gcx*ssy;
        ERP[i][1] = Gcx*ssz;
        ERP[i][2] = Gcy*ssz;
        ERP[i][3] = 0.5*(xssxww - yssyww);
        ERP[i][4] = SQR3I*(Gcz*ssz
                           - 0.5*(xssxww + yssyww));
    }
}

/*****************************************************************************/
/*  Primitive 3-center OVP evaluation [ ABC ] for type [ DPS ]           */
/*  np      ; number of first primitive pair data            */
/*  Zp      ; (Alpha+Beta) from orbital exponents Alpha, Beta    */
/*  px, py, pz  ; (Alpha*ax * Beta*bx) / (Alpha+Beta) , ...          */
/*  pax, pay, paz   ; ax, ay, az                         */
/*  pbx, pby, pbz   ; bx, by, bz                         */
/*  HP      ; Zeta**(-3/2) * exp(-(Alpha*Beta/Zeta) * (A-B)**2   */
/*  Gamma       ; orbital exponent Gamma for last GTO            */
/*  cx, cy, cz  ; coordinates of last GTO C              */
/*  cc      ; coefficient of last GTO C              */
/*  ERP     ; output primitive OVP                   */
/*              ERP[0,1,2,...,np-1][0] for [ DxyPxS ]        */
/*              ERP[0,1,2,...,np-1][1] for [ DxyPyS ]        */
/*              ERP[0,1,2,...,np-1][2] for [ DxyPzS ]        */
/*              ERP[0,1,2,...,np-1][3] for [ DxzPxS ]        */
/*              ERP[0,1,2,...,np-1][4] for [ DxzPyS ]        */
/*              ERP[0,1,2,...,np-1][5] for [ DxzPzS ]        */
/*              ERP[0,1,2,...,np-1][6] for [ DyzPxS ]        */
/*              ERP[0,1,2,...,np-1][7] for [ DyzPyS ]        */
/*              ERP[0,1,2,...,np-1][8] for [ DyzPzS ]        */
/*              ERP[0,1,2,...,np-1][9] for [ Dxx-yyPxS ]     */
/*                  .............                */
/*              ERP[0,1,2,...,np-1][12] for [ D3zz-rrPxS ]   */
/*                  .............                */
/*****************************************************************************/
void DfOverlap::ovpDPS(const int np,
                       const double* px, const double* py, const double* pz,
                       const double* pax, const double* pay, const double* paz,
                       const double* pbx, const double* pby, const double* pbz,
                       const double Gamma,
                       const double cx, const double cy, const double cz,
                       const double cc,
                       const double* Zp, const double* HP, double** ERP)
{
    const double hc  = M_PI*sqrt(M_PI)*cc;
    for (int i = 0; i < np; ++i) {
        const double RZG = 1.0 / (Zp[i] + Gamma);
        const double Roo = Zp[i]*Gamma*RZG;
        const double pcx = px[i] - cx;
        const double pcy = py[i] - cy;
        const double pcz = pz[i] - cz;
        const double t   = Roo * (pcx*pcx + pcy*pcy + pcz*pcz);
        const double sss = hc * HP[i] * RZG * sqrt(RZG) * exp(-t);

        const double Gx  = (Zp[i]*px[i] + Gamma*cx) * RZG;
        const double Gy  = (Zp[i]*py[i] + Gamma*cy) * RZG;
        const double Gz  = (Zp[i]*pz[i] + Gamma*cz) * RZG;
        const double Gax = Gx - pax[i];
        const double Gay = Gy - pay[i];
        const double Gaz = Gz - paz[i];
        const double Gbx = Gx - pbx[i];
        const double Gby = Gy - pby[i];
        const double Gbz = Gz - pbz[i];
        // [ PSS ]
        const double RZGP5 = 0.5*RZG;
        const double xss   = Gax*sss;
        const double yss   = Gay*sss;
        const double zss   = Gaz*sss;
        // [ DSS ]  xy, xz, yz, xx-yy, 3zz-rr
        const double xxssww = Gax*xss;
        const double yyssww = Gay*yss;

        const double ass = Gax*yss;
        const double bss = Gax*zss;
        const double css = Gay*zss;
        const double dss = 0.5*(xxssww - yyssww);
        const double ess = SQR3I*(Gaz*zss
                                  - 0.5*(xxssww + yyssww));
        // [ DPS ]
        const double Rxss = RZGP5*xss;
        const double Ryss = RZGP5*yss;
        const double Rzss = RZGP5*zss;

        ERP[i][ 0] = Gbx*ass + Ryss;
        ERP[i][ 1] = Gby*ass + Rxss;
        ERP[i][ 2] = Gbz*ass;
        ERP[i][ 3] = Gbx*bss + Rzss;
        ERP[i][ 4] = Gby*bss;
        ERP[i][ 5] = Gbz*bss + Rxss;
        ERP[i][ 6] = Gbx*css;
        ERP[i][ 7] = Gby*css + Rzss;
        ERP[i][ 8] = Gbz*css + Ryss;
        ERP[i][ 9] = Gbx*dss + Rxss;
        ERP[i][10] = Gby*dss - Ryss;
        ERP[i][11] = Gbz*dss;
        ERP[i][12] = Gbx*ess - SQR3I *Rxss;
        ERP[i][13] = Gby*ess - SQR3I *Ryss;
        ERP[i][14] = Gbz*ess + SQR3I2*Rzss;
    }
}

/*****************************************************************************/
/*  Primitive 3-center OVP evaluation [ ABC ] for type [ DSP ]           */
/*  np      ; number of first primitive pair data            */
/*  Zp      ; (Alpha+Beta) from orbital exponents Alpha, Beta    */
/*  px, py, pz  ; (Alpha*ax * Beta*bx) / (Alpha+Beta) , ...          */
/*  pax, pay, paz   ; ax, ay, az                         */
/*  HP      ; Zeta**(-3/2) * exp(-(Alpha*Beta/Zeta) * (A-B)**2   */
/*  Gamma       ; orbital exponent Gamma for last GTO            */
/*  cx, cy, cz  ; coordinates of last GTO C              */
/*  cc      ; coefficient of last GTO C              */
/*  ERP     ; output primitive OVP                   */
/*              ERP[0,1,2,...,np-1][0] for [ DxySPx ]        */
/*              ERP[0,1,2,...,np-1][1] for [ DxySPy ]        */
/*              ERP[0,1,2,...,np-1][2] for [ DxySPz ]        */
/*              ERP[0,1,2,...,np-1][3] for [ DxzSPx ]        */
/*              ERP[0,1,2,...,np-1][4] for [ DxzSPy ]        */
/*              ERP[0,1,2,...,np-1][5] for [ DxzSPz ]        */
/*              ERP[0,1,2,...,np-1][6] for [ DyzSPx ]        */
/*              ERP[0,1,2,...,np-1][7] for [ DyzSPy ]        */
/*              ERP[0,1,2,...,np-1][8] for [ DyzSPz ]        */
/*              ERP[0,1,2,...,np-1][9] for [ Dxx-yySPx ]     */
/*                  .............                */
/*              ERP[0,1,2,...,np-1][12] for [ D3zz-rrSPx ]   */
/*                  .............                */
/*****************************************************************************/
void DfOverlap::ovpDSP(const int np,
                       const double* px, const double* py, const double* pz,
                       const double* pax, const double* pay, const double* paz,
                       const double Gamma,
                       const double cx, const double cy, const double cz,
                       const double cc,
                       const double* Zp, const double* HP, double** ERP)
{
    const double hc  = M_PI*sqrt(M_PI)*cc;
    for (int i = 0; i < np; ++i) {
        const double RZG = 1.0 / (Zp[i] + Gamma);
        const double Roo = Zp[i]*Gamma*RZG;
        const double pcx = px[i] - cx;
        const double pcy = py[i] - cy;
        const double pcz = pz[i] - cz;
        const double t   = Roo * (pcx*pcx + pcy*pcy + pcz*pcz);
        const double sss = hc * HP[i] * RZG * sqrt(RZG) * exp(-t);

        const double Gx  = (Zp[i]*px[i] + Gamma*cx) * RZG;
        const double Gy  = (Zp[i]*py[i] + Gamma*cy) * RZG;
        const double Gz  = (Zp[i]*pz[i] + Gamma*cz) * RZG;
        const double Gax = Gx - pax[i];
        const double Gay = Gy - pay[i];
        const double Gaz = Gz - paz[i];
        const double Gcx = Gx - cx;
        const double Gcy = Gy - cy;
        const double Gcz = Gz - cz;
        // [ PSS ]
        const double RZGP5 = 0.5*RZG;
        const double xss   = Gax*sss;
        const double yss   = Gay*sss;
        const double zss   = Gaz*sss;
        // [ DSS ]  xy, xz, yz, xx-yy, 3zz-rr
        const double xxssww = Gax*xss;
        const double yyssww = Gay*yss;

        const double ass = Gax*yss;
        const double bss = Gax*zss;
        const double css = Gay*zss;
        const double dss = 0.5*(xxssww - yyssww);
        const double ess = SQR3I*(Gaz*zss
                                  - 0.5*(xxssww + yyssww));
        // [ DSP ]
        const double Rxss = RZGP5*xss;
        const double Ryss = RZGP5*yss;
        const double Rzss = RZGP5*zss;

        ERP[i][ 0] = Gcx*ass + Ryss;
        ERP[i][ 1] = Gcy*ass + Rxss;
        ERP[i][ 2] = Gcz*ass;
        ERP[i][ 3] = Gcx*bss + Rzss;
        ERP[i][ 4] = Gcy*bss;
        ERP[i][ 5] = Gcz*bss + Rxss;
        ERP[i][ 6] = Gcx*css;
        ERP[i][ 7] = Gcy*css + Rzss;
        ERP[i][ 8] = Gcz*css + Ryss;
        ERP[i][ 9] = Gcx*dss + Rxss;
        ERP[i][10] = Gcy*dss - Ryss;
        ERP[i][11] = Gcz*dss;
        ERP[i][12] = Gcx*ess - SQR3I *Rxss;
        ERP[i][13] = Gcy*ess - SQR3I *Ryss;
        ERP[i][14] = Gcz*ess + SQR3I2*Rzss;
    }
}

/*****************************************************************************/
/*                                       */
/*  Primitive 3-center OVP evaluation [ ABC ] for type [ PSD ]               */
/*                                       */
/*  np      ; number of first primitive pair data            */
/*  Zp      ; (Alpha+Beta) from orbital exponents Alpha, Beta    */
/*  px, py, pz  ; (Alpha*ax * Beta*bx) / (Alpha+Beta) , ...          */
/*  pax, pay, paz   ; ax, ay, az                         */
/*  HP      ; Zeta**(-3/2) * exp(-(Alpha*Beta/Zeta) * (A-B)**2   */
/*  Gamma       ; orbital exponent Gamma for last GTO            */
/*  cx, cy, cz  ; coordinates of last GTO C              */
/*  cc      ; coefficient of last GTO C              */
/*  ERP     ; output primitive OVP                   */
/*              ERP[0,1,2,...,np-1][0] for [ PxSDxy ]        */
/*              ERP[0,1,2,...,np-1][1] for [ PxSDxz ]        */
/*              ERP[0,1,2,...,np-1][2] for [ PxSDyz ]        */
/*              ERP[0,1,2,...,np-1][3] for [ PxSDxx-yy ]     */
/*              ERP[0,1,2,...,np-1][4] for [ PxSD3zz-rr ]    */
/*              ERP[0,1,2,...,np-1][5] for [ PySDxy ]        */
/*              ERP[0,1,2,...,np-1][6] for [ PySDxz ]        */
/*              ERP[0,1,2,...,np-1][7] for [ PySDyz ]        */
/*              ERP[0,1,2,...,np-1][8] for [ PySDxx-yy]      */
/*              ERP[0,1,2,...,np-1][9] for [ PySD3zz-rr ]    */
/*                  .............                */
/*                                       */
/***************************************************************************/
void DfOverlap::ovpPSD(const int np,
                       const double* px, const double* py, const double* pz,
                       const double* pax, const double* pay, const double* paz,
                       const double Gamma,
                       const double cx, const double cy, const double cz,
                       const double cc,
                       const double* Zp, const double* HP, double** ERP)
{
    const double hc = M_PI*sqrt(M_PI)*cc;
    for (int i = 0; i < np; ++i) {
        const double RZG = 1.0 / (Zp[i] + Gamma);
        const double Roo = Zp[i]*Gamma*RZG;
        const double pcx = px[i] - cx;
        const double pcy = py[i] - cy;
        const double pcz = pz[i] - cz;
        const double t   = Roo * (pcx*pcx + pcy*pcy + pcz*pcz);
        const double sss = hc * HP[i] * RZG * sqrt(RZG) * exp(-t);

        const double Gx  = (Zp[i]*px[i] + Gamma*cx) * RZG;
        const double Gy  = (Zp[i]*py[i] + Gamma*cy) * RZG;
        const double Gz  = (Zp[i]*pz[i] + Gamma*cz) * RZG;
        const double Gax = Gx - pax[i];
        const double Gay = Gy - pay[i];
        const double Gaz = Gz - paz[i];
        const double Gcx = Gx - cx;
        const double Gcy = Gy - cy;
        const double Gcz = Gz - cz;
        // [ SSP ]
        const double RZGP5 = 0.5*RZG;
        const double ssx   = Gcx*sss;
        const double ssy   = Gcy*sss;
        const double ssz   = Gcz*sss;
        // [ SSD ]  xy, xz, yz, xx-yy, 3zz-rr
        const double xssxww = Gcx*ssx;
        const double yssyww = Gcy*ssy;

        const double ssa = Gcx*ssy;
        const double ssb = Gcx*ssz;
        const double ssc = Gcy*ssz;
        const double ssd = 0.5*(xssxww - yssyww);
        const double sse = SQR3I*(Gcz*ssz
                                  - 0.5*(xssxww + yssyww));
        // [ PSD ]
        const double Rssx = RZGP5*ssx;
        const double Rssy = RZGP5*ssy;
        const double Rssz = RZGP5*ssz;

        ERP[i][ 0] = Gax*ssa + Rssy;
        ERP[i][ 1] = Gax*ssb + Rssz;
        ERP[i][ 2] = Gax*ssc;
        ERP[i][ 3] = Gax*ssd + Rssx;
        ERP[i][ 4] = Gax*sse - SQR3I *Rssx;
        ERP[i][ 5] = Gay*ssa + Rssx;
        ERP[i][ 6] = Gay*ssb;
        ERP[i][ 7] = Gay*ssc + Rssz;
        ERP[i][ 8] = Gay*ssd - Rssy;
        ERP[i][ 9] = Gay*sse - SQR3I *Rssy;
        ERP[i][10] = Gaz*ssa;
        ERP[i][11] = Gaz*ssb + Rssx;
        ERP[i][12] = Gaz*ssc + Rssy;
        ERP[i][13] = Gaz*ssd;
        ERP[i][14] = Gaz*sse + SQR3I2*Rssz;
    }
}


/*****************************************************************************/
/*                                       */
/*  Primitive 3-center OVP evaluation [ ABC ] for type [ DPP ]           */
/*                                       */
/*  np      ; number of first primitive pair data            */
/*  Zp      ; (Alpha+Beta) from orbital exponents Alpha, Beta    */
/*  px, py, pz  ; (Alpha*ax * Beta*bx) / (Alpha+Beta) , ...          */
/*  pax, pay, paz   ; ax, ay, az                         */
/*  pbx, pby, pbz   ; bx, by, bz                         */
/*  HP      ; Zeta**(-3/2) * exp(-(Alpha*Beta/Zeta) * (A-B)**2   */
/*  Gamma       ; orbital exponent Gamma for last GTO            */
/*  cx, cy, cz  ; coordinates of last GTO C              */
/*  cc      ; coefficient of last GTO C              */
/*  ERP     ; output primitive OVP                   */
/*              ERP[0,1,2,...,np-1][0] for [ DxyPxPx ]       */
/*              ERP[0,1,2,...,np-1][1] for [ DxyPxPy ]       */
/*              ERP[0,1,2,...,np-1][2] for [ DxyPxPz ]       */
/*              ERP[0,1,2,...,np-1][3] for [ DxyPyPx ]       */
/*              ERP[0,1,2,...,np-1][4] for [ DxyPyPy ]       */
/*              ERP[0,1,2,...,np-1][5] for [ DxyPyPz ]       */
/*              ERP[0,1,2,...,np-1][6] for [ DxyPzPx ]       */
/*              ERP[0,1,2,...,np-1][7] for [ DxyPzPy ]       */
/*              ERP[0,1,2,...,np-1][8] for [ DxyPzPz ]       */
/*                  .............                */
/*                                       */
/*****************************************************************************/
void DfOverlap::ovpDPP(const int np,
                       const double* px, const double* py, const double* pz,
                       const double* pax, const double* pay, const double* paz,
                       const double* pbx, const double* pby, const double* pbz,
                       const double Gamma,
                       const double cx, const double cy, const double cz,
                       const double cc,
                       const double* Zp, const double* HP, double** ERP)
{
    const double hc = M_PI*sqrt(M_PI)*cc;
    for (int i = 0; i < np; ++i) {
        const double RZG = 1.0 / (Zp[i] + Gamma);
        const double Roo = Zp[i]*Gamma*RZG;
        const double pcx = px[i] - cx;
        const double pcy = py[i] - cy;
        const double pcz = pz[i] - cz;
        const double t   = Roo * (pcx*pcx + pcy*pcy + pcz*pcz);
        const double sss = hc * HP[i] * RZG * sqrt(RZG) * exp(-t);

        const double Gx  = (Zp[i]*px[i] + Gamma*cx) * RZG;
        const double Gy  = (Zp[i]*py[i] + Gamma*cy) * RZG;
        const double Gz  = (Zp[i]*pz[i] + Gamma*cz) * RZG;
        const double Gax = Gx - pax[i];
        const double Gay = Gy - pay[i];
        const double Gaz = Gz - paz[i];
        const double Gbx = Gx - pbx[i];
        const double Gby = Gy - pby[i];
        const double Gbz = Gz - pbz[i];
        const double Gcx = Gx - cx;
        const double Gcy = Gy - cy;
        const double Gcz = Gz - cz;
        // [ PSS ]
        const double RZGP5 = 0.5*RZG;
        const double xss   = Gax*sss;
        const double yss   = Gay*sss;
        const double zss   = Gaz*sss;
        // [ DSS ]  xy, xz, yz, xx-yy, 3zz-rr
        const double xxssww = Gax*xss;
        const double yyssww = Gay*yss;

        const double ass = Gax*yss;
        const double bss = Gax*zss;
        const double css = Gay*zss;
        const double dss = 0.5*(xxssww - yyssww);
        const double ess = SQR3I*(Gaz*zss
                                  - 0.5*(xxssww + yyssww));
        // [ PPS ]
        const double Rsss = RZGP5*sss;
        const double xxs  = Gbx*xss + Rsss;
        const double xys  = Gby*xss;
        const double xzs  = Gbz*xss;
        const double yxs  = Gbx*yss;
        const double yys  = Gby*yss + Rsss;
        const double yzs  = Gbz*yss;
        const double zxs  = Gbx*zss;
        const double zys  = Gby*zss;
        const double zzs  = Gbz*zss + Rsss;
        // [ DPS ]
        const double Rxss = RZGP5*xss;
        const double Ryss = RZGP5*yss;
        const double Rzss = RZGP5*zss;
        const double axs  = Gbx*ass + Ryss;
        const double ays  = Gby*ass + Rxss;
        const double azs  = Gbz*ass;
        const double bxs  = Gbx*bss + Rzss;
        const double bys  = Gby*bss;
        const double bzs  = Gbz*bss + Rxss;
        const double cxs  = Gbx*css;
        const double cys  = Gby*css + Rzss;
        const double czs  = Gbz*css + Ryss;
        const double dxs  = Gbx*dss + Rxss;
        const double dys  = Gby*dss - Ryss;
        const double dzs  = Gbz*dss;
        const double exs  = Gbx*ess - SQR3I *Rxss;
        const double eys  = Gby*ess - SQR3I *Ryss;
        const double ezs  = Gbz*ess + SQR3I2*Rzss;
        // [ DPP ]
        const double Rass = RZGP5*ass;
        const double Rbss = RZGP5*bss;
        const double Rcss = RZGP5*css;
        const double Rdss = RZGP5*dss;
        const double Ress = RZGP5*ess;
        const double Rxxs = RZGP5*xxs;
        const double Rxys = RZGP5*xys;
        const double Rxzs = RZGP5*xzs;
        const double Ryxs = RZGP5*yxs;
        const double Ryys = RZGP5*yys;
        const double Ryzs = RZGP5*yzs;
        const double Rzxs = RZGP5*zxs;
        const double Rzys = RZGP5*zys;
        const double Rzzs = RZGP5*zzs;

        // a = xy
        ERP[i][ 0] = Gcx*axs + Ryxs + Rass;
        ERP[i][ 1] = Gcy*axs + Rxxs;
        ERP[i][ 2] = Gcz*axs;
        ERP[i][ 3] = Gcx*ays + Ryys;
        ERP[i][ 4] = Gcy*ays + Rxys + Rass;
        ERP[i][ 5] = Gcz*ays;
        ERP[i][ 6] = Gcx*azs + Ryzs;
        ERP[i][ 7] = Gcy*azs + Rxzs;
        ERP[i][ 8] = Gcz*azs        + Rass;
        // b = xz
        ERP[i][ 9] = Gcx*bxs + Rzxs + Rbss;
        ERP[i][10] = Gcy*bxs;
        ERP[i][11] = Gcz*bxs + Rxxs;
        ERP[i][12] = Gcx*bys + Rzys;
        ERP[i][13] = Gcy*bys        + Rbss;
        ERP[i][14] = Gcz*bys + Rxys;
        ERP[i][15] = Gcx*bzs + Rzzs;
        ERP[i][16] = Gcy*bzs;
        ERP[i][17] = Gcz*bzs + Rxzs + Rbss;
        // c = yz
        ERP[i][18] = Gcx*cxs        + Rcss;
        ERP[i][19] = Gcy*cxs + Rzxs;
        ERP[i][20] = Gcz*cxs + Ryxs;
        ERP[i][21] = Gcx*cys;
        ERP[i][22] = Gcy*cys + Rzys + Rcss;
        ERP[i][23] = Gcz*cys + Ryys;
        ERP[i][24] = Gcx*czs;
        ERP[i][25] = Gcy*czs + Rzzs;
        ERP[i][26] = Gcz*czs + Ryzs + Rcss;
        // d = xx-yy
        ERP[i][27] = Gcx*dxs + Rxxs + Rdss;
        ERP[i][28] = Gcy*dxs - Ryxs;
        ERP[i][29] = Gcz*dxs;
        ERP[i][30] = Gcx*dys + Rxys;
        ERP[i][31] = Gcy*dys - Ryys + Rdss;
        ERP[i][32] = Gcz*dys;
        ERP[i][33] = Gcx*dzs + Rxzs;
        ERP[i][34] = Gcy*dzs - Ryzs;
        ERP[i][35] = Gcz*dzs        + Rdss;
        // e = 3zz-rr
        ERP[i][36] = Gcx*exs - SQR3I *Rxxs + Ress;
        ERP[i][37] = Gcy*exs - SQR3I *Ryxs;
        ERP[i][38] = Gcz*exs + SQR3I2*Rzxs;
        ERP[i][39] = Gcx*eys - SQR3I *Rxys;
        ERP[i][40] = Gcy*eys - SQR3I *Ryys + Ress;
        ERP[i][41] = Gcz*eys + SQR3I2*Rzys;
        ERP[i][42] = Gcx*ezs - SQR3I *Rxzs;
        ERP[i][43] = Gcy*ezs - SQR3I *Ryzs;
        ERP[i][44] = Gcz*ezs + SQR3I2*Rzzs + Ress;
    }
}


/*****************************************************************************/
/*                                       */
/*  Primitive 3-center OVP evaluation [ ABC ] for type [ PPD ]           */
/*                                       */
/*  np      ; number of first primitive pair data            */
/*  Zp      ; (Alpha+Beta) from orbital exponents Alpha, Beta    */
/*  px, py, pz  ; (Alpha*ax * Beta*bx) / (Alpha+Beta) , ...          */
/*  pax, pay, paz   ; ax, ay, az                         */
/*  pbx, pby, pbz   ; bx, by, bz                         */
/*  HP      ; Zeta**(-3/2) * exp(-(Alpha*Beta/Zeta) * (A-B)**2   */
/*  Gamma       ; orbital exponent Gamma for last GTO            */
/*  cx, cy, cz  ; coordinates of last GTO C              */
/*  cc      ; coefficient of last GTO C              */
/*  ERP     ; output primitive OVP                   */
/*              ERP[0,1,2,...,np-1][0] for [ PxPxDxy ]       */
/*              ERP[0,1,2,...,np-1][1] for [ PxPxDxz ]       */
/*              ERP[0,1,2,...,np-1][2] for [ PxPxDyz ]       */
/*              ERP[0,1,2,...,np-1][3] for [ PxPxDxx-yy ]    */
/*              ERP[0,1,2,...,np-1][4] for [ PxPxD3zz-rr ]   */
/*              ERP[0,1,2,...,np-1][5] for [ PxPyDxy ]       */
/*              ERP[0,1,2,...,np-1][6] for [ PxPyDxz ]       */
/*              ERP[0,1,2,...,np-1][7] for [ PxPyDyz ]       */
/*              ERP[0,1,2,...,np-1][8] for [ PxPyDxx-yy ]    */
/*              ERP[0,1,2,...,np-1][8] for [ PxPyD3zz-rr ]   */
/*                  .............                */
/*                                       */
/*****************************************************************************/
void DfOverlap::ovpPPD(const int np,
                       const double* px, const double* py, const double* pz,
                       const double* pax, const double* pay, const double* paz,
                       const double* pbx, const double* pby, const double* pbz,
                       const double Gamma,
                       const double cx, const double cy, const double cz,
                       const double cc,
                       const double* Zp, const double* HP, double** ERP)
{
    const double hc  = M_PI*sqrt(M_PI)*cc;
    for (int i = 0; i < np; ++i) {
        const double RZG = 1.0 / (Zp[i] + Gamma);
        const double Roo = Zp[i]*Gamma*RZG;
        const double pcx = px[i] - cx;
        const double pcy = py[i] - cy;
        const double pcz = pz[i] - cz;
        const double t   = Roo * (pcx*pcx + pcy*pcy + pcz*pcz);
        const double sss = hc * HP[i] * RZG * sqrt(RZG) * exp(-t);

        const double Gx  = (Zp[i]*px[i] + Gamma*cx) * RZG;
        const double Gy  = (Zp[i]*py[i] + Gamma*cy) * RZG;
        const double Gz  = (Zp[i]*pz[i] + Gamma*cz) * RZG;
        const double Gax = Gx - pax[i];
        const double Gay = Gy - pay[i];
        const double Gaz = Gz - paz[i];
        const double Gbx = Gx - pbx[i];
        const double Gby = Gy - pby[i];
        const double Gbz = Gz - pbz[i];
        const double Gcx = Gx - cx;
        const double Gcy = Gy - cy;
        const double Gcz = Gz - cz;
        // [ SSP ]
        const double RZGP5 = 0.5*RZG;
        const double ssx   = Gcx*sss;
        const double ssy   = Gcy*sss;
        const double ssz   = Gcz*sss;
        // [ PSP ]
        const double RGF = RZGP5*sss;
        const double xsx = Gax*ssx + RGF;
        const double xsy = Gax*ssy;
        const double xsz = Gax*ssz;
        const double ysx = Gay*ssx;
        const double ysy = Gay*ssy + RGF;
        const double ysz = Gay*ssz;
        const double zsx = Gaz*ssx;
        const double zsy = Gaz*ssy;
        const double zsz = Gaz*ssz + RGF;
        // [ SSD ]  xy, xz, yz, xx-yy, 3zz-rr
        const double xssxww = Gcx*ssx;
        const double yssyww = Gcy*ssy;

        const double ssa = Gcx*ssy;
        const double ssb = Gcx*ssz;
        const double ssc = Gcy*ssz;
        const double ssd = 0.5*(xssxww - yssyww);
        const double sse = SQR3I*(Gcz*ssz
                                  - 0.5*(xssxww + yssyww));
        // [ PSD ]
        const double Rssx = RZGP5*ssx;
        const double Rssy = RZGP5*ssy;
        const double Rssz = RZGP5*ssz;
        const double xsa = Gax*ssa + Rssy;
        const double xsb = Gax*ssb + Rssz;
        const double xsc = Gax*ssc;
        const double xsd = Gax*ssd + Rssx;
        const double xse = Gax*sse - SQR3I *Rssx;
        const double ysa = Gay*ssa + Rssx;
        const double ysb = Gay*ssb;
        const double ysc = Gay*ssc + Rssz;
        const double ysd = Gay*ssd - Rssy;
        const double yse = Gay*sse - SQR3I *Rssy;
        const double zsa = Gaz*ssa;
        const double zsb = Gaz*ssb + Rssx;
        const double zsc = Gaz*ssc + Rssy;
        const double zsd = Gaz*ssd;
        const double zse = Gaz*sse + SQR3I2*Rssz;
        // [ PPD ]
        const double Rssa = RZGP5*ssa;
        const double Rssb = RZGP5*ssb;
        const double Rssc = RZGP5*ssc;
        const double Rssd = RZGP5*ssd;
        const double Rsse = RZGP5*sse;
        const double Rxsx = RZGP5*xsx;
        const double Rxsy = RZGP5*xsy;
        const double Rxsz = RZGP5*xsz;
        const double Rysx = RZGP5*ysx;
        const double Rysy = RZGP5*ysy;
        const double Rysz = RZGP5*ysz;
        const double Rzsx = RZGP5*zsx;
        const double Rzsy = RZGP5*zsy;
        const double Rzsz = RZGP5*zsz;

        ERP[i][ 0] = Gbx*xsa + Rssa + Rxsy;
        ERP[i][ 1] = Gbx*xsb + Rssb + Rxsz;
        ERP[i][ 2] = Gbx*xsc + Rssc;
        ERP[i][ 3] = Gbx*xsd + Rssd + Rxsx;
        ERP[i][ 4] = Gbx*xse + Rsse - SQR3I *Rxsx;
        ERP[i][ 5] = Gby*xsa        + Rxsx;
        ERP[i][ 6] = Gby*xsb;
        ERP[i][ 7] = Gby*xsc        + Rxsz;
        ERP[i][ 8] = Gby*xsd        - Rxsy;
        ERP[i][ 9] = Gby*xse        - SQR3I *Rxsy;
        ERP[i][10] = Gbz*xsa;
        ERP[i][11] = Gbz*xsb        + Rxsx;
        ERP[i][12] = Gbz*xsc        + Rxsy;
        ERP[i][13] = Gbz*xsd;
        ERP[i][14] = Gbz*xse        + SQR3I2*Rxsz;

        ERP[i][15] = Gbx*ysa        + Rysy;
        ERP[i][16] = Gbx*ysb        + Rysz;
        ERP[i][17] = Gbx*ysc;
        ERP[i][18] = Gbx*ysd        + Rysx;
        ERP[i][19] = Gbx*yse        - SQR3I *Rysx;
        ERP[i][20] = Gby*ysa + Rssa + Rysx;
        ERP[i][21] = Gby*ysb + Rssb;
        ERP[i][22] = Gby*ysc + Rssc + Rysz;
        ERP[i][23] = Gby*ysd + Rssd - Rysy;
        ERP[i][24] = Gby*yse + Rsse - SQR3I *Rysy;
        ERP[i][25] = Gbz*ysa;
        ERP[i][26] = Gbz*ysb        + Rysx;
        ERP[i][27] = Gbz*ysc        + Rysy;
        ERP[i][28] = Gbz*ysd;
        ERP[i][29] = Gbz*yse        + SQR3I2*Rysz;

        ERP[i][30] = Gbx*zsa        + Rzsy;
        ERP[i][31] = Gbx*zsb        + Rzsz;
        ERP[i][32] = Gbx*zsc;
        ERP[i][33] = Gbx*zsd        + Rzsx;
        ERP[i][34] = Gbx*zse        - SQR3I *Rzsx;
        ERP[i][35] = Gby*zsa        + Rzsx;
        ERP[i][36] = Gby*zsb;
        ERP[i][37] = Gby*zsc        + Rzsz;
        ERP[i][38] = Gby*zsd        - Rzsy;
        ERP[i][39] = Gby*zse        - SQR3I *Rzsy;
        ERP[i][40] = Gbz*zsa + Rssa;
        ERP[i][41] = Gbz*zsb + Rssb + Rzsx;
        ERP[i][42] = Gbz*zsc + Rssc + Rzsy;
        ERP[i][43] = Gbz*zsd + Rssd;
        ERP[i][44] = Gbz*zse + Rsse + SQR3I2*Rzsz;
    }
}


/*****************************************************************************/
/*                                       */
/*  Primitive 3-center OVP evaluation [ ABC ] for type [ DDS ]           */
/*                                       */
/*  np      ; number of first primitive pair data            */
/*  Zp      ; (Alpha+Beta) from orbital exponents Alpha, Beta    */
/*  px, py, pz  ; (Alpha*ax * Beta*bx) / (Alpha+Beta) , ...          */
/*  pax, pay, paz   ; ax, ay, az                         */
/*  pbx, pby, pbz   ; bx, by, bz                         */
/*  HP      ; Zeta**(-3/2) * exp(-(Alpha*Beta/Zeta) * (A-B)**2   */
/*  Gamma       ; orbital exponent Gamma for last GTO            */
/*  cx, cy, cz  ; coordinates of last GTO C              */
/*  cc      ; coefficient of last GTO C              */
/*  ERP     ; output primitive OVP                   */
/*              ERP[0,1,2,...,np-1][0] for [ DxyDxyS ]       */
/*              ERP[0,1,2,...,np-1][1] for [ DxyDxzS ]       */
/*              ERP[0,1,2,...,np-1][2] for [ DxyDyzS ]       */
/*              ERP[0,1,2,...,np-1][3] for [ DxyDxx-yyS ]    */
/*              ERP[0,1,2,...,np-1][4] for [ DxyD3zz-rrS ]   */
/*              ERP[0,1,2,...,np-1][5] for [ DxzDxyS ]       */
/*              ERP[0,1,2,...,np-1][6] for [ DxzDxzS ]       */
/*              ERP[0,1,2,...,np-1][7] for [ DxzDyzS ]       */
/*              ERP[0,1,2,...,np-1][8] for [ DxzDxx-yyS ]    */
/*              ERP[0,1,2,...,np-1][9] for [ DxzD3zz-rrS ]   */
/*                  .............                */
/*                                       */
/*****************************************************************************/
void DfOverlap::ovpDDS(const int np,
                       const double* px, const double* py, const double* pz,
                       const double* pax, const double* pay, const double* paz,
                       const double* pbx, const double* pby, const double* pbz,
                       const double Gamma,
                       const double cx, const double cy, const double cz,
                       const double cc,
                       const double* Zp, const double* HP, double** ERP)
{
    const double hc = M_PI*sqrt(M_PI)*cc;
    for (int i = 0; i < np; ++i) {
        const double RZG = 1.0 / (Zp[i] + Gamma);
        const double Roo = Zp[i]*Gamma*RZG;
        const double pcx = px[i] - cx;
        const double pcy = py[i] - cy;
        const double pcz = pz[i] - cz;
        const double t   = Roo * (pcx*pcx + pcy*pcy + pcz*pcz);
        const double sss = hc * HP[i] * RZG * sqrt(RZG) * exp(-t);

        const double Gx  = (Zp[i]*px[i] + Gamma*cx) * RZG;
        const double Gy  = (Zp[i]*py[i] + Gamma*cy) * RZG;
        const double Gz  = (Zp[i]*pz[i] + Gamma*cz) * RZG;
        const double Gax = Gx - pax[i];
        const double Gay = Gy - pay[i];
        const double Gaz = Gz - paz[i];
        const double Gbx = Gx - pbx[i];
        const double Gby = Gy - pby[i];
        const double Gbz = Gz - pbz[i];
        // [ PSS ]
        const double RZGP5 = 0.5*RZG;
        const double xss   = Gax*sss;
        const double yss   = Gay*sss;
        const double zss   = Gaz*sss;
        // [ PPS ]
        const double Rsss = RZGP5*sss;
        const double xxs  = Gbx*xss + Rsss;
        const double xys  = Gby*xss;
        const double xzs  = Gbz*xss;
        const double yxs  = Gbx*yss;
        const double yys  = Gby*yss + Rsss;
        const double yzs  = Gbz*yss;
        const double zxs  = Gbx*zss;
        const double zys  = Gby*zss;
        const double zzs  = Gbz*zss + Rsss;
        // [ DSS ]  xy, xz, yz, xx-yy, 3zz-rr
        const double xxssww = Gax*xss;
        const double yyssww = Gay*yss;

        const double ass = Gax*yss;
        const double bss = Gax*zss;
        const double css = Gay*zss;
        const double dss = 0.5*(xxssww - yyssww);
        const double ess = SQR3I*(Gaz*zss
                                  - 0.5*(xxssww + yyssww));
        // [ DPS ]
        const double Rxss = RZGP5*xss;
        const double Ryss = RZGP5*yss;
        const double Rzss = RZGP5*zss;
        const double axs = Gbx*ass + Ryss;
        const double ays = Gby*ass + Rxss;
        const double azs = Gbz*ass;
        const double bxs = Gbx*bss + Rzss;
        const double bys = Gby*bss;
        const double bzs = Gbz*bss + Rxss;
        const double cxs = Gbx*css;
        const double cys = Gby*css + Rzss;
        const double czs = Gbz*css + Ryss;
        const double dxs = Gbx*dss + Rxss;
        const double dys = Gby*dss - Ryss;
        const double dzs = Gbz*dss;
        const double exs = Gbx*ess - SQR3I *Rxss;
        const double eys = Gby*ess - SQR3I *Ryss;
        const double ezs = Gbz*ess + SQR3I2*Rzss;
        // [ DDS ]
        const double Rxxs  = RZGP5*xxs;
        const double Rxys  = RZGP5*xys;
        const double Rxzs  = RZGP5*xzs;
        const double Ryxs  = RZGP5*yxs;
        const double Ryys  = RZGP5*yys;
        const double Ryzs  = RZGP5*yzs;
        const double Rzxs  = RZGP5*zxs;
        const double Rzys  = RZGP5*zys;
        const double Rzzs  = RZGP5*zzs;
        const double axxsw = Gbx*axs + Ryxs;
        const double ayysw = Gby*ays + Rxys;
        const double azzsw = Gbz*azs;
        const double bxxsw = Gbx*bxs + Rzxs;
        const double byysw = Gby*bys;
        const double bzzsw = Gbz*bzs + Rxzs;
        const double cxxsw = Gbx*cxs;
        const double cyysw = Gby*cys + Rzys;
        const double czzsw = Gbz*czs + Ryzs;
        const double dxxsw = Gbx*dxs + Rxxs;
        const double dyysw = Gby*dys - Ryys;
        const double dzzsw = Gbz*dzs;
        const double exxsw = Gbx*exs - SQR3I *Rxxs;
        const double eyysw = Gby*eys - SQR3I *Ryys;
        const double ezzsw = Gbz*ezs + SQR3I2*Rzzs;

        ERP[i][ 0] = Gbx*ays + Ryys;
        ERP[i][ 1] = Gbx*azs + Ryzs;
        ERP[i][ 2] = Gby*azs + Rxzs;
        ERP[i][ 3] = 0.5*(axxsw - ayysw);
        ERP[i][ 4] = SQR3I*(azzsw - 0.5*(axxsw + ayysw));
        ERP[i][ 5] = Gbx*bys + Rzys;
        ERP[i][ 6] = Gbx*bzs + Rzzs;
        ERP[i][ 7] = Gby*bzs;
        ERP[i][ 8] = 0.5*(bxxsw - byysw);
        ERP[i][ 9] = SQR3I*(bzzsw - 0.5*(bxxsw + byysw));
        ERP[i][10] = Gbx*cys;
        ERP[i][11] = Gbx*czs;
        ERP[i][12] = Gby*czs + Rzzs;
        ERP[i][13] = 0.5*(cxxsw - cyysw);
        ERP[i][14] = SQR3I*(czzsw - 0.5*(cxxsw + cyysw));
        ERP[i][15] = Gbx*dys + Rxys;
        ERP[i][16] = Gbx*dzs + Rxzs;
        ERP[i][17] = Gby*dzs - Ryzs;
        ERP[i][18] = 0.5*(dxxsw - dyysw);
        ERP[i][19] = SQR3I*(dzzsw - 0.5*(dxxsw + dyysw));
        ERP[i][20] = Gbx*eys - SQR3I*Rxys;
        ERP[i][21] = Gbx*ezs - SQR3I*Rxzs;
        ERP[i][22] = Gby*ezs - SQR3I*Ryzs;
        ERP[i][23] = 0.5*(exxsw - eyysw);
        ERP[i][24] = SQR3I*(ezzsw - 0.5*(exxsw + eyysw));
    }
}


/*****************************************************************************/
/*                                       */
/*  Primitive 3-center OVP evaluation [ ABC ] for type [ DSD ]           */
/*                                       */
/*  np      ; number of first primitive pair data            */
/*  Zp      ; (Alpha+Beta) from orbital exponents Alpha, Beta    */
/*  px, py, pz  ; (Alpha*ax * Beta*bx) / (Alpha+Beta) , ...          */
/*  pax, pay, paz   ; ax, ay, az                         */
/*  HP      ; Zeta**(-3/2) * exp(-(Alpha*Beta/Zeta) * (A-B)**2   */
/*  Gamma       ; orbital exponent Gamma for last GTO            */
/*  cx, cy, cz  ; coordinates of last GTO C              */
/*  cc      ; coefficient of last GTO C              */
/*  ERP     ; output primitive OVP                   */
/*              ERP[0,1,2,...,np-1][0] for [ DxySDxy ]       */
/*              ERP[0,1,2,...,np-1][1] for [ DxySDxz ]       */
/*              ERP[0,1,2,...,np-1][2] for [ DxySDyz ]       */
/*              ERP[0,1,2,...,np-1][3] for [ DxySDxx-yy ]    */
/*              ERP[0,1,2,...,np-1][4] for [ DxySD3zz-rr ]   */
/*              ERP[0,1,2,...,np-1][5] for [ DxzSDxy ]       */
/*              ERP[0,1,2,...,np-1][6] for [ DxzSDxz ]       */
/*              ERP[0,1,2,...,np-1][7] for [ DxzSDyz ]       */
/*              ERP[0,1,2,...,np-1][8] for [ DxzSDxx-yy ]    */
/*              ERP[0,1,2,...,np-1][9] for [ DxzSD3zz-rr ]   */
/*                  .............                */
/*                                       */
/*****************************************************************************/
void DfOverlap::ovpDSD(const int np,
                       const double* px, const double* py, const double* pz,
                       const double* pax, const double* pay, const double* paz,
                       const double Gamma,
                       const double cx, const double cy, const double cz,
                       const double cc,
                       const double* Zp, const double* HP, double** ERP)
{
    const double hc = M_PI*sqrt(M_PI)*cc;
    for (int i = 0; i < np; ++i) {
        const double RZG = 1.0 / (Zp[i] + Gamma);
        const double Roo = Zp[i]*Gamma*RZG;
        const double pcx = px[i] - cx;
        const double pcy = py[i] - cy;
        const double pcz = pz[i] - cz;
        const double t   = Roo * (pcx*pcx + pcy*pcy + pcz*pcz);
        const double sss = hc * HP[i] * RZG * sqrt(RZG) * exp(-t);

        const double Gx  = (Zp[i]*px[i] + Gamma*cx) * RZG;
        const double Gy  = (Zp[i]*py[i] + Gamma*cy) * RZG;
        const double Gz  = (Zp[i]*pz[i] + Gamma*cz) * RZG;
        const double Gax = Gx - pax[i];
        const double Gay = Gy - pay[i];
        const double Gaz = Gz - paz[i];
        const double Gcx = Gx - cx;
        const double Gcy = Gy - cy;
        const double Gcz = Gz - cz;
        // [ PSS ]
        const double RZGP5 = 0.5*RZG;
        const double xss   = Gax*sss;
        const double yss   = Gay*sss;
        const double zss   = Gaz*sss;
        // [ PSP ]
        const double RGF = RZGP5*sss;
        const double xsx = Gcx*xss + RGF;
        const double xsy = Gcy*xss;
        const double xsz = Gcz*xss;
        const double ysx = Gcx*yss;
        const double ysy = Gcy*yss + RGF;
        const double ysz = Gcz*yss;
        const double zsx = Gcx*zss;
        const double zsy = Gcy*zss;
        const double zsz = Gcz*zss + RGF;
        // [ DSS ]  xy, xz, yz, xx-yy, 3zz-rr
        const double xxssww = Gax*xss;
        const double yyssww = Gay*yss;
        const double ass = Gax*yss;
        const double bss = Gax*zss;
        const double css = Gay*zss;
        const double dss = 0.5*(xxssww - yyssww);
        const double ess = SQR3I*(Gaz*zss
                                  - 0.5*(xxssww + yyssww));
        // [ DSP ]
        const double Rxss = RZGP5*xss;
        const double Ryss = RZGP5*yss;
        const double Rzss = RZGP5*zss;
        const double asx = Gcx*ass + Ryss;
        const double asy = Gcy*ass + Rxss;
        const double asz = Gcz*ass;
        const double bsx = Gcx*bss + Rzss;
        const double bsy = Gcy*bss;
        const double bsz = Gcz*bss + Rxss;
        const double csx = Gcx*css;
        const double csy = Gcy*css + Rzss;
        const double csz = Gcz*css + Ryss;
        const double dsx = Gcx*dss + Rxss;
        const double dsy = Gcy*dss - Ryss;
        const double dsz = Gcz*dss;
        const double esx = Gcx*ess - SQR3I *Rxss;
        const double esy = Gcy*ess - SQR3I *Ryss;
        const double esz = Gcz*ess + SQR3I2*Rzss;
        // [ DSD ]
        const double Rxsx  = RZGP5*xsx;
        const double Rxsy  = RZGP5*xsy;
        const double Rxsz  = RZGP5*xsz;
        const double Rysx  = RZGP5*ysx;
        const double Rysy  = RZGP5*ysy;
        const double Rysz  = RZGP5*ysz;
        const double Rzsx  = RZGP5*zsx;
        const double Rzsy  = RZGP5*zsy;
        const double Rzsz  = RZGP5*zsz;
        const double asxxw = Gcx*asx + Rysx;
        const double asyyw = Gcy*asy + Rxsy;
        const double bsxxw = Gcx*bsx + Rzsx;
        const double bsyyw = Gcy*bsy;
        const double csxxw = Gcx*csx;
        const double csyyw = Gcy*csy + Rzsy;
        const double dsxxw = Gcx*dsx + Rxsx;
        const double dsyyw = Gcy*dsy - Rysy;
        const double esxxw = Gcx*esx - SQR3I *Rxsx;
        const double esyyw = Gcy*esy - SQR3I *Rysy;

        ERP[i][ 0] = Gcx*asy + Rysy;
        ERP[i][ 1] = Gcx*asz + Rysz;
        ERP[i][ 2] = Gcy*asz + Rxsz;
        ERP[i][ 3] = 0.5*(asxxw - asyyw);
        ERP[i][ 4] = SQR3I*(Gcz*asz
                            - 0.5*(asxxw + asyyw));
        ERP[i][ 5] = Gcx*bsy + Rzsy;
        ERP[i][ 6] = Gcx*bsz + Rzsz;
        ERP[i][ 7] = Gcy*bsz;
        ERP[i][ 8] = 0.5*(bsxxw - bsyyw);
        ERP[i][ 9] = SQR3I*(Gcz*bsz + Rxsz
                            - 0.5*(bsxxw + bsyyw));
        ERP[i][10] = Gcx*csy;
        ERP[i][11] = Gcx*csz;
        ERP[i][12] = Gcy*csz + Rzsz;
        ERP[i][13] = 0.5*(csxxw - csyyw);
        ERP[i][14] = SQR3I*(Gcz*csz + Rysz
                            - 0.5*(csxxw + csyyw));
        ERP[i][15] = Gcx*dsy + Rxsy;
        ERP[i][16] = Gcx*dsz + Rxsz;
        ERP[i][17] = Gcy*dsz - Rysz;
        ERP[i][18] = 0.5*(dsxxw - dsyyw);
        ERP[i][19] = SQR3I*(Gcz*dsz
                            - 0.5*(dsxxw + dsyyw));
        ERP[i][20] = Gcx*esy - SQR3I*Rxsy;
        ERP[i][21] = Gcx*esz - SQR3I*Rxsz;
        ERP[i][22] = Gcy*esz - SQR3I*Rysz;
        ERP[i][23] = 0.5*(esxxw - esyyw);
        ERP[i][24] = SQR3I*(Gcz*esz + SQR3I2*Rzsz
                            - 0.5*(esxxw + esyyw));
    }
}


/*****************************************************************************/
/*                                       */
/*  Primitive 3-center OVP evaluation [ ABC ] for type [ DDP ]           */
/*                                       */
/*  np      ; number of first primitive pair data            */
/*  Zp      ; (Alpha+Beta) from orbital exponents Alpha, Beta    */
/*  px, py, pz  ; (Alpha*ax * Beta*bx) / (Alpha+Beta) , ...          */
/*  pax, pay, paz   ; ax, ay, az                         */
/*  pbx, pby, pbz   ; bx, by, bz                         */
/*  HP      ; Zeta**(-3/2) * exp(-(Alpha*Beta/Zeta) * (A-B)**2   */
/*  Gamma       ; orbital exponent Gamma for last GTO            */
/*  cx, cy, cz  ; coordinates of last GTO C              */
/*  cc      ; coefficient of last GTO C              */
/*  ERP     ; output primitive OVP                   */
/*              ERP[0,1,2,...,np-1][0] for [ DxyDxyPx ]      */
/*              ERP[0,1,2,...,np-1][1] for [ DxyDxyPy ]      */
/*              ERP[0,1,2,...,np-1][2] for [ DxyDxyPz ]      */
/*              ERP[0,1,2,...,np-1][3] for [ DxyDxzPx ]      */
/*              ERP[0,1,2,...,np-1][4] for [ DxyDxzPy ]      */
/*              ERP[0,1,2,...,np-1][5] for [ DxyDxzPz ]      */
/*              ERP[0,1,2,...,np-1][6] for [ DxyDyzPx ]      */
/*              ERP[0,1,2,...,np-1][7] for [ DxyDyzPy ]      */
/*              ERP[0,1,2,...,np-1][8] for [ DxyDyzPz ]      */
/*              ERP[0,1,2,...,np-1][9] for [ DxyDxx-yyPx ]   */
/*                  .............                */
/*                                       */
/*****************************************************************************/
void DfOverlap::ovpDDP(const int np,
                       const double* px, const double* py, const double* pz,
                       const double* pax, const double* pay, const double* paz,
                       const double* pbx, const double* pby, const double* pbz,
                       const double Gamma,
                       const double cx, const double cy, const double cz,
                       const double cc,
                       const double* Zp, const double* HP, double** ERP)
{
    const double hc = M_PI*sqrt(M_PI)*cc;
    for (int i = 0; i < np; ++i) {
        const double RZG = 1.0 / (Zp[i] + Gamma);
        const double Roo = Zp[i]*Gamma*RZG;
        const double pcx = px[i] - cx;
        const double pcy = py[i] - cy;
        const double pcz = pz[i] - cz;
        const double t   = Roo * (pcx*pcx + pcy*pcy + pcz*pcz);
        const double sss = hc * HP[i] * RZG * sqrt(RZG) * exp(-t);

        const double Gx  = (Zp[i]*px[i] + Gamma*cx) * RZG;
        const double Gy  = (Zp[i]*py[i] + Gamma*cy) * RZG;
        const double Gz  = (Zp[i]*pz[i] + Gamma*cz) * RZG;
        const double Gax = Gx - pax[i];
        const double Gay = Gy - pay[i];
        const double Gaz = Gz - paz[i];
        const double Gbx = Gx - pbx[i];
        const double Gby = Gy - pby[i];
        const double Gbz = Gz - pbz[i];
        const double Gcx = Gx - cx;
        const double Gcy = Gy - cy;
        const double Gcz = Gz - cz;
        const double abx = pax[i] - pbx[i];
        const double aby = pay[i] - pby[i];
        const double abz = paz[i] - pbz[i];
        // [ PSS ]
        const double RZGP5 = 0.5*RZG;
        const double xss   = Gax*sss;
        const double yss   = Gay*sss;
        const double zss   = Gaz*sss;
        // [ PPS ]
        const double Rsss = RZGP5*sss;
        const double xxs  = Gbx*xss + Rsss;
        const double xys  = Gby*xss;
        const double xzs  = Gbz*xss;
        const double yxs  = Gbx*yss;
        const double yys  = Gby*yss + Rsss;
        const double yzs  = Gbz*yss;
        const double zxs  = Gbx*zss;
        const double zys  = Gby*zss;
        const double zzs  = Gbz*zss + Rsss;
        // [ DSS ]  xy, xz, yz, xx-yy, 3zz-rr
        const double xxssww = Gax*xss;
        const double yyssww = Gay*yss;
        const double ass = Gax*yss;
        const double bss = Gax*zss;
        const double css = Gay*zss;
        const double dss = 0.5*(xxssww - yyssww);
        const double ess = SQR3I*(Gaz*zss
                                  - 0.5*(xxssww + yyssww));
        // [ DSS ] -2
        const double x2ss = xxssww  + Rsss;
        const double y2ss = yyssww  + Rsss;
        const double z2ss = Gaz*zss + Rsss;
        // [ DPS ]
        const double Rxss = RZGP5*xss;
        const double Ryss = RZGP5*yss;
        const double Rzss = RZGP5*zss;
        const double axs = Gbx*ass + Ryss;
        const double ays = Gby*ass + Rxss;
        const double azs = Gbz*ass;
        const double bxs = Gbx*bss + Rzss;
        const double bys = Gby*bss;
        const double bzs = Gbz*bss + Rxss;
        const double cxs = Gbx*css;
        const double cys = Gby*css + Rzss;
        const double czs = Gbz*css + Ryss;
        const double dxs = Gbx*dss + Rxss;
        const double dys = Gby*dss - Ryss;
        const double dzs = Gbz*dss;
        const double exs = Gbx*ess - SQR3I *Rxss;
        const double eys = Gby*ess - SQR3I *Ryss;
        const double ezs = Gbz*ess + SQR3I2*Rzss;
        // [ DPS ] -2 (Dxx)xs
        const double x2xs = Gbx*x2ss + Rxss + Rxss;
        const double y2ys = Gby*y2ss + Ryss + Ryss;
        const double z2zs = Gbz*z2ss + Rzss + Rzss;
        // [ PDS ]
        const double xas = axs + aby*xxs;
        const double xbs = bxs + abz*xxs;
        const double xcs = azs + aby*xzs;
        const double xds = 0.5*(x2xs + abx*xxs - ays - aby*xys);
        const double xes = SQR3I*(bzs + abz*xzs
                                  -0.5*(x2xs + abx*xxs + ays + aby*xys));
        const double yas = ays + abx*yys;
        const double ybs = cxs + abz*yxs;
        const double ycs = cys + abz*yys;
        const double yds = 0.5*(axs + abx*yxs - y2ys - aby*yys);
        const double yes = SQR3I*(czs + abz*yzs
                                  -0.5*(axs + abx*yxs + y2ys + aby*yys));
        const double zas = bys + abx*zys;
        const double zbs = bzs + abx*zzs;
        const double zcs = czs + aby*zzs;
        const double zds = 0.5*(bxs + abx*zxs - cys - aby*zys);
        const double zes = SQR3I*(z2zs + abz*zzs
                                  -0.5*(bxs + abx*zxs + cys + aby*zys));
        // [ DPP ]
        //const double Rass = RZGP5*ass;
        //const double Rbss = RZGP5*bss;
        //const double Rcss = RZGP5*css;
        //const double Rdss = RZGP5*dss;
        //const double Ress = RZGP5*ess;
        const double Rxxs = RZGP5*xxs;
        const double Rxys = RZGP5*xys;
        const double Rxzs = RZGP5*xzs;
        const double Ryxs = RZGP5*yxs;
        const double Ryys = RZGP5*yys;
        const double Ryzs = RZGP5*yzs;
        const double Rzxs = RZGP5*zxs;
        const double Rzys = RZGP5*zys;
        const double Rzzs = RZGP5*zzs;
        // a = xy
        //const double axx = Gcx*axs + Ryxs + Rass;
        //const double axy = Gcy*axs + Rxxs;
        //const double axz = Gcz*axs;
        //const double ayx = Gcx*ays + Ryys;
        //const double ayy = Gcy*ays + Rxys + Rass;
        //const double ayz = Gcz*ays;
        //const double azx = Gcx*azs + Ryzs;
        //const double azy = Gcy*azs + Rxzs;
        //const double azz = Gcz*azs        + Rass;
        // b = xz
        //const double bxx = Gcx*bxs + Rzxs + Rbss;
        //const double bxy = Gcy*bxs;
        //const double bxz = Gcz*bxs + Rxxs;
        //const double byx = Gcx*bys + Rzys;
        //const double byy = Gcy*bys        + Rbss;
        //const double byz = Gcz*bys + Rxys;
        //const double bzx = Gcx*bzs + Rzzs;
        //const double bzy = Gcy*bzs;
        //const double bzz = Gcz*bzs + Rxzs + Rbss;
        // c = yz
        //const double cxx = Gcx*cxs        + Rcss;
        //const double cxy = Gcy*cxs + Rzxs;
        //const double cxz = Gcz*cxs + Ryxs;
        //const double cyx = Gcx*cys;
        //const double cyy = Gcy*cys + Rzys + Rcss;
        //const double cyz = Gcz*cys + Ryys;
        //const double czx = Gcx*czs;
        //const double czy = Gcy*czs + Rzzs;
        //const double czz = Gcz*czs + Ryzs + Rcss;
        // d = xx-yy
        //const double dxx = Gcx*dxs + Rxxs + Rdss;
        //const double dxy = Gcy*dxs - Ryxs;
        //const double dxz = Gcz*dxs;
        //const double dyx = Gcx*dys + Rxys;
        //const double dyy = Gcy*dys - Ryys + Rdss;
        //const double dyz = Gcz*dys;
        //const double dzx = Gcx*dzs + Rxzs;
        //const double dzy = Gcy*dzs - Ryzs;
        //const double dzz = Gcz*dzs        + Rdss;
        // e = 3zz-rr
        //const double exx = Gcx*exs - SQR3I *Rxxs + Ress;
        //const double exy = Gcy*exs - SQR3I *Ryxs;
        //const double exz = Gcz*exs + SQR3I2*Rzxs;
        //const double eyx = Gcx*eys - SQR3I *Rxys;
        //const double eyy = Gcy*eys - SQR3I *Ryys + Ress;
        //const double eyz = Gcz*eys + SQR3I2*Rzys;
        //const double ezx = Gcx*ezs - SQR3I *Rxzs;
        //const double ezy = Gcy*ezs - SQR3I *Ryzs;
        //const double ezz = Gcz*ezs + SQR3I2*Rzzs + Ress;
        // [ DDS ]
        //const double Rxxs  = RZGP5*xxs;
        //const double Rxys  = RZGP5*xys;
        //const double Rxzs  = RZGP5*xzs;
        //const double Ryxs  = RZGP5*yxs;
        //const double Ryys  = RZGP5*yys;
        //const double Ryzs  = RZGP5*yzs;
        //const double Rzxs  = RZGP5*zxs;
        //const double Rzys  = RZGP5*zys;
        //const double Rzzs  = RZGP5*zzs;
        const double axxsw = Gbx*axs + Ryxs;
        const double ayysw = Gby*ays + Rxys;
        const double azzsw = Gbz*azs;
        const double bxxsw = Gbx*bxs + Rzxs;
        const double byysw = Gby*bys;
        const double bzzsw = Gbz*bzs + Rxzs;
        const double cxxsw = Gbx*cxs;
        const double cyysw = Gby*cys + Rzys;
        const double czzsw = Gbz*czs + Ryzs;
        const double dxxsw = Gbx*dxs + Rxxs;
        const double dyysw = Gby*dys - Ryys;
        const double dzzsw = Gbz*dzs;
        const double exxsw = Gbx*exs - SQR3I *Rxxs;
        const double eyysw = Gby*eys - SQR3I *Ryys;
        const double ezzsw = Gbz*ezs + SQR3I2*Rzzs;
        const double aas = Gbx*ays + Ryys;
        const double abs = Gbx*azs + Ryzs;
        const double acs = Gby*azs + Rxzs;
        const double ads = 0.5*(axxsw - ayysw);
        const double aes = SQR3I*(azzsw - 0.5*(axxsw + ayysw));
        const double bas = Gbx*bys + Rzys;
        const double bbs = Gbx*bzs + Rzzs;
        const double bcs = Gby*bzs;
        const double bds = 0.5*(bxxsw - byysw);
        const double bes = SQR3I*(bzzsw - 0.5*(bxxsw + byysw));
        const double cas = Gbx*cys;
        const double cbs = Gbx*czs;
        const double ccs = Gby*czs + Rzzs;
        const double cds = 0.5*(cxxsw - cyysw);
        const double ces = SQR3I*(czzsw - 0.5*(cxxsw + cyysw));
        const double das = Gbx*dys + Rxys;
        const double dbs = Gbx*dzs + Rxzs;
        const double dcs = Gby*dzs - Ryzs;
        const double dds = 0.5*(dxxsw - dyysw);
        const double des = SQR3I*(dzzsw - 0.5*(dxxsw + dyysw));
        const double eas = Gbx*eys - SQR3I*Rxys;
        const double ebs = Gbx*ezs - SQR3I*Rxzs;
        const double ecs = Gby*ezs - SQR3I*Ryzs;
        const double eds = 0.5*(exxsw - eyysw);
        const double ees = SQR3I*(ezzsw - 0.5*(exxsw + eyysw));
        // [ DDP ]
        const double Raxs = RZGP5*axs;
        const double Rays = RZGP5*ays;
        const double Razs = RZGP5*azs;
        const double Rbxs = RZGP5*bxs;
        const double Rbys = RZGP5*bys;
        const double Rbzs = RZGP5*bzs;
        const double Rcxs = RZGP5*cxs;
        const double Rcys = RZGP5*cys;
        const double Rczs = RZGP5*czs;
        const double Rdxs = RZGP5*dxs;
        const double Rdys = RZGP5*dys;
        const double Rdzs = RZGP5*dzs;
        const double Rexs = RZGP5*exs;
        const double Reys = RZGP5*eys;
        const double Rezs = RZGP5*ezs;
        const double Rxas = RZGP5*xas;
        const double Rxbs = RZGP5*xbs;
        const double Rxcs = RZGP5*xcs;
        const double Rxds = RZGP5*xds;
        const double Rxes = RZGP5*xes;
        const double Ryas = RZGP5*yas;
        const double Rybs = RZGP5*ybs;
        const double Rycs = RZGP5*ycs;
        const double Ryds = RZGP5*yds;
        const double Ryes = RZGP5*yes;
        const double Rzas = RZGP5*zas;
        const double Rzbs = RZGP5*zbs;
        const double Rzcs = RZGP5*zcs;
        const double Rzds = RZGP5*zds;
        const double Rzes = RZGP5*zes;

        ERP[i][ 0] = Gcx*aas + Ryas + Rays;
        ERP[i][ 1] = Gcy*aas + Rxas + Raxs;
        ERP[i][ 2] = Gcz*aas;
        ERP[i][ 3] = Gcx*abs + Rybs + Razs;
        ERP[i][ 4] = Gcy*abs + Rxbs;
        ERP[i][ 5] = Gcz*abs        + Raxs;
        ERP[i][ 6] = Gcx*acs + Rycs;
        ERP[i][ 7] = Gcy*acs + Rxcs + Razs;
        ERP[i][ 8] = Gcz*acs        + Rays;
        ERP[i][ 9] = Gcx*ads + Ryds + Raxs;
        ERP[i][10] = Gcy*ads + Rxds - Rays;
        ERP[i][11] = Gcz*ads;
        ERP[i][12] = Gcx*aes + Ryes - SQR3I *Raxs;
        ERP[i][13] = Gcy*aes + Rxes - SQR3I *Rays;
        ERP[i][14] = Gcz*aes        + SQR3I2*Razs;

        ERP[i][15] = Gcx*bas + Rzas + Rbys;
        ERP[i][16] = Gcy*bas        + Rbxs;
        ERP[i][17] = Gcz*bas + Rxas;
        ERP[i][18] = Gcx*bbs + Rzbs + Rbzs;
        ERP[i][19] = Gcy*bbs;
        ERP[i][20] = Gcz*bbs + Rxbs + Rbxs;
        ERP[i][21] = Gcx*bcs + Rzcs;
        ERP[i][22] = Gcy*bcs        + Rbzs;
        ERP[i][23] = Gcz*bcs + Rxcs + Rbys;
        ERP[i][24] = Gcx*bds + Rzds + Rbxs;
        ERP[i][25] = Gcy*bds        - Rbys;
        ERP[i][26] = Gcz*bds + Rxds;
        ERP[i][27] = Gcx*bes + Rzes - SQR3I *Rbxs;
        ERP[i][28] = Gcy*bes        - SQR3I *Rbys;
        ERP[i][29] = Gcz*bes + Rxes + SQR3I2*Rbzs;

        ERP[i][30] = Gcx*cas        + Rcys;
        ERP[i][31] = Gcy*cas + Rzas + Rcxs;
        ERP[i][32] = Gcz*cas + Ryas;
        ERP[i][33] = Gcx*cbs        + Rczs;
        ERP[i][34] = Gcy*cbs + Rzbs;
        ERP[i][35] = Gcz*cbs + Rybs + Rcxs;
        ERP[i][36] = Gcx*ccs;
        ERP[i][37] = Gcy*ccs + Rzcs + Rczs;
        ERP[i][38] = Gcz*ccs + Rycs + Rcys;
        ERP[i][39] = Gcx*cds        + Rcxs;
        ERP[i][40] = Gcy*cds + Rzds - Rcys;
        ERP[i][41] = Gcz*cds + Ryds;
        ERP[i][42] = Gcx*ces        - SQR3I *Rcxs;
        ERP[i][43] = Gcy*ces + Rzes - SQR3I *Rcys;
        ERP[i][44] = Gcz*ces + Ryes + SQR3I2*Rczs;

        ERP[i][45] = Gcx*das + Rxas + Rdys;
        ERP[i][46] = Gcy*das - Ryas + Rdxs;
        ERP[i][47] = Gcz*das;
        ERP[i][48] = Gcx*dbs + Rxbs + Rdzs;
        ERP[i][49] = Gcy*dbs - Rybs;
        ERP[i][50] = Gcz*dbs        + Rdxs;
        ERP[i][51] = Gcx*dcs + Rxcs;
        ERP[i][52] = Gcy*dcs - Rycs + Rdzs;
        ERP[i][53] = Gcz*dcs        + Rdys;
        ERP[i][54] = Gcx*dds + Rxds + Rdxs;
        ERP[i][55] = Gcy*dds - Ryds - Rdys;
        ERP[i][56] = Gcz*dds;
        ERP[i][57] = Gcx*des + Rxes - SQR3I *Rdxs;
        ERP[i][58] = Gcy*des - Ryes - SQR3I *Rdys;
        ERP[i][59] = Gcz*des        + SQR3I2*Rdzs;

        ERP[i][60] = Gcx*eas - SQR3I *Rxas + Reys;
        ERP[i][61] = Gcy*eas - SQR3I *Ryas + Rexs;
        ERP[i][62] = Gcz*eas + SQR3I2*Rzas;
        ERP[i][63] = Gcx*ebs - SQR3I *Rxbs + Rezs;
        ERP[i][64] = Gcy*ebs - SQR3I *Rybs;
        ERP[i][65] = Gcz*ebs + SQR3I2*Rzbs + Rexs;
        ERP[i][66] = Gcx*ecs - SQR3I *Rxcs;
        ERP[i][67] = Gcy*ecs - SQR3I *Rycs + Rezs;
        ERP[i][68] = Gcz*ecs + SQR3I2*Rzcs + Reys;
        ERP[i][69] = Gcx*eds - SQR3I *Rxds + Rexs;
        ERP[i][70] = Gcy*eds - SQR3I *Ryds - Reys;
        ERP[i][71] = Gcz*eds + SQR3I2*Rzds;
        ERP[i][72] = Gcx*ees - SQR3I *(Rxes + Rexs);
        ERP[i][73] = Gcy*ees - SQR3I *(Ryes + Reys);
        ERP[i][74] = Gcz*ees + SQR3I2*(Rzes + Rezs);
    }
}


/*****************************************************************************/
/*                                       */
/*  Primitive 3-center OVP evaluation [ ABC ] for type [ DPD ]           */
/*                                       */
/*  np      ; number of first primitive pair data            */
/*  Zp      ; (Alpha+Beta) from orbital exponents Alpha, Beta    */
/*  px, py, pz  ; (Alpha*ax * Beta*bx) / (Alpha+Beta) , ...          */
/*  pax, pay, paz   ; ax, ay, az                         */
/*  pbx, pby, pbz   ; bx, by, bz                         */
/*  HP      ; Zeta**(-3/2) * exp(-(Alpha*Beta/Zeta) * (A-B)**2   */
/*  Gamma       ; orbital exponent Gamma for last GTO            */
/*  cx, cy, cz  ; coordinates of last GTO C              */
/*  cc      ; coefficient of last GTO C              */
/*  ERP     ; output primitive OVP                   */
/*          ERP[0,1,2,...,np-1][0] for [ DxyPxDxy ]          */
/*          ERP[0,1,2,...,np-1][1] for [ DxyPxDxz ]          */
/*          ERP[0,1,2,...,np-1][2] for [ DxyPxDyz ]          */
/*          ERP[0,1,2,...,np-1][3] for [ DxyPxDxx-yy ]       */
/*          ERP[0,1,2,...,np-1][4] for [ DxyPxD3zz-rr ]      */
/*          ERP[0,1,2,...,np-1][5] for [ DxyPyDxy ]          */
/*          ERP[0,1,2,...,np-1][6] for [ DxyPyDxz ]          */
/*          ERP[0,1,2,...,np-1][7] for [ DxyPyDyz ]          */
/*          ERP[0,1,2,...,np-1][8] for [ DxyPyDxx-yy ]           */
/*          ERP[0,1,2,...,np-1][9] for [ DxyPyD3zz-rr ]      */
/*                  .............                */
/*                                       */
/*****************************************************************************/
void DfOverlap::ovpDPD(const int np,
                       const double* px, const double* py, const double* pz,
                       const double* pax, const double* pay, const double* paz,
                       const double* pbx, const double* pby, const double* pbz,
                       const double Gamma,
                       const double cx, const double cy, const double cz,
                       const double cc,
                       const double* Zp, const double* HP, double** ERP)
{
    const double hc = M_PI*sqrt(M_PI)*cc;
    for (int i = 0; i < np; ++i) {
        const double RZG = 1.0 / (Zp[i] + Gamma);
        const double Roo = Zp[i]*Gamma*RZG;
        const double pcx = px[i] - cx;
        const double pcy = py[i] - cy;
        const double pcz = pz[i] - cz;
        const double t   = Roo * (pcx*pcx + pcy*pcy + pcz*pcz);
        const double sss = hc * HP[i] * RZG * sqrt(RZG) * exp(-t);

        const double Gx  = (Zp[i]*px[i] + Gamma*cx) * RZG;
        const double Gy  = (Zp[i]*py[i] + Gamma*cy) * RZG;
        const double Gz  = (Zp[i]*pz[i] + Gamma*cz) * RZG;
        const double Gax = Gx - pax[i];
        const double Gay = Gy - pay[i];
        const double Gaz = Gz - paz[i];
        const double Gbx = Gx - pbx[i];
        const double Gby = Gy - pby[i];
        const double Gbz = Gz - pbz[i];
        const double Gcx = Gx - cx;
        const double Gcy = Gy - cy;
        const double Gcz = Gz - cz;
        // [ PSS ]
        const double RZGP5 = 0.5*RZG;
        const double xss   = Gax*sss;
        const double yss   = Gay*sss;
        const double zss   = Gaz*sss;
        // [ SSP ]
        const double ssx   = Gcx*sss;
        const double ssy   = Gcy*sss;
        const double ssz   = Gcz*sss;
        // [ PSP ]
        const double RGF = RZGP5*sss;
        const double xsx = Gcx*xss + RGF;
        const double xsy = Gcy*xss;
        const double xsz = Gcz*xss;
        const double ysx = Gcx*yss;
        const double ysy = Gcy*yss + RGF;
        const double ysz = Gcz*yss;
        const double zsx = Gcx*zss;
        const double zsy = Gcy*zss;
        const double zsz = Gcz*zss + RGF;
        // [ DSS ]  xy, xz, yz, xx-yy, 3zz-rr
        const double xxssww = Gax*xss;
        const double yyssww = Gay*yss;
        const double ass = Gax*yss;
        const double bss = Gax*zss;
        const double css = Gay*zss;
        const double dss = 0.5*(xxssww - yyssww);
        const double ess = SQR3I*(Gaz*zss
                                  - 0.5*(xxssww + yyssww));
        // [ DSP ]
        const double Rxss = RZGP5*xss;
        const double Ryss = RZGP5*yss;
        const double Rzss = RZGP5*zss;
        const double asx = Gcx*ass + Ryss;
        const double asy = Gcy*ass + Rxss;
        const double asz = Gcz*ass;
        const double bsx = Gcx*bss + Rzss;
        const double bsy = Gcy*bss;
        const double bsz = Gcz*bss + Rxss;
        const double csx = Gcx*css;
        const double csy = Gcy*css + Rzss;
        const double csz = Gcz*css + Ryss;
        const double dsx = Gcx*dss + Rxss;
        const double dsy = Gcy*dss - Ryss;
        const double dsz = Gcz*dss;
        const double esx = Gcx*ess - SQR3I *Rxss;
        const double esy = Gcy*ess - SQR3I *Ryss;
        const double esz = Gcz*ess + SQR3I2*Rzss;
        // [ PSD ]
        const double xsx2w = Gcx*xsx + RZGP5*ssx;
        const double ysx2w = Gcx*ysx;
        const double zsx2w = Gcx*zsx;
        const double xsy2w = Gcy*xsy;
        const double ysy2w = Gcy*ysy + RZGP5*ssy;
        const double zsy2w = Gcy*zsy;
        const double xsz2w = Gcz*xsz;
        const double ysz2w = Gcz*ysz;
        const double zsz2w = Gcz*zsz + RZGP5*ssz;
        const double xsa   = Gcy*xsx;
        const double ysa   = Gcx*ysy;
        const double zsa   = Gcx*zsy;
        const double xsb   = Gcz*xsx;
        const double ysb   = Gcx*ysz;
        const double zsb   = Gcx*zsz;
        const double xsc   = Gcy*xsz;
        const double ysc   = Gcz*ysy;
        const double zsc   = Gcy*zsz;
        const double xsd   = 0.5*(xsx2w - xsy2w);
        const double ysd   = 0.5*(ysx2w - ysy2w);
        const double zsd   = 0.5*(zsx2w - zsy2w);
        const double xse   = SQR3I*(xsz2w - 0.5*(xsx2w + xsy2w));
        const double yse   = SQR3I*(ysz2w - 0.5*(ysx2w + ysy2w));
        const double zse   = SQR3I*(zsz2w - 0.5*(zsx2w + zsy2w));
        // [ DSD ]
        const double Rxsx  = RZGP5*xsx;
        const double Rxsy  = RZGP5*xsy;
        const double Rxsz  = RZGP5*xsz;
        const double Rysx  = RZGP5*ysx;
        const double Rysy  = RZGP5*ysy;
        const double Rysz  = RZGP5*ysz;
        const double Rzsx  = RZGP5*zsx;
        const double Rzsy  = RZGP5*zsy;
        const double Rzsz  = RZGP5*zsz;
        const double asxxw = Gcx*asx + Rysx;
        const double asyyw = Gcy*asy + Rxsy;
        const double bsxxw = Gcx*bsx + Rzsx;
        const double bsyyw = Gcy*bsy;
        const double csxxw = Gcx*csx;
        const double csyyw = Gcy*csy + Rzsy;
        const double dsxxw = Gcx*dsx + Rxsx;
        const double dsyyw = Gcy*dsy - Rysy;
        const double esxxw = Gcx*esx - SQR3I *Rxsx;
        const double esyyw = Gcy*esy - SQR3I *Rysy;

        const double asa = Gcx*asy + Rysy;
        const double asb = Gcx*asz + Rysz;
        const double asc = Gcy*asz + Rxsz;
        const double asd = 0.5*(asxxw - asyyw);
        const double ase = SQR3I*(Gcz*asz
                                  - 0.5*(asxxw + asyyw));
        const double bsa = Gcx*bsy + Rzsy;
        const double bsb = Gcx*bsz + Rzsz;
        const double bsc = Gcy*bsz;
        const double bsd = 0.5*(bsxxw - bsyyw);
        const double bse = SQR3I*(Gcz*bsz + Rxsz
                                  - 0.5*(bsxxw + bsyyw));
        const double csa = Gcx*csy;
        const double csb = Gcx*csz;
        const double csc = Gcy*csz + Rzsz;
        const double csd = 0.5*(csxxw - csyyw);
        const double cse = SQR3I*(Gcz*csz + Rysz
                                  - 0.5*(csxxw + csyyw));
        const double dsa = Gcx*dsy + Rxsy;
        const double dsb = Gcx*dsz + Rxsz;
        const double dsc = Gcy*dsz - Rysz;
        const double dsd = 0.5*(dsxxw - dsyyw);
        const double dse = SQR3I*(Gcz*dsz
                                  - 0.5*(dsxxw + dsyyw));
        const double esa = Gcx*esy - SQR3I*Rxsy;
        const double esb = Gcx*esz - SQR3I*Rxsz;
        const double esc = Gcy*esz - SQR3I*Rysz;
        const double esd = 0.5*(esxxw - esyyw);
        const double ese = SQR3I*(Gcz*esz + SQR3I2*Rzsz
                                  - 0.5*(esxxw + esyyw));
        // [ DPD ]
        const double Rxsa = RZGP5*xsa;
        const double Rysa = RZGP5*ysa;
        const double Rzsa = RZGP5*zsa;
        const double Rxsb = RZGP5*xsb;
        const double Rysb = RZGP5*ysb;
        const double Rzsb = RZGP5*zsb;
        const double Rxsc = RZGP5*xsc;
        const double Rysc = RZGP5*ysc;
        const double Rzsc = RZGP5*zsc;
        const double Rxsd = RZGP5*xsd;
        const double Rysd = RZGP5*ysd;
        const double Rzsd = RZGP5*zsd;
        const double Rxse = RZGP5*xse;
        const double Ryse = RZGP5*yse;
        const double Rzse = RZGP5*zse;
        const double Rasx = RZGP5*asx;
        const double Rasy = RZGP5*asy;
        const double Rasz = RZGP5*asz;
        const double Rbsx = RZGP5*bsx;
        const double Rbsy = RZGP5*bsy;
        const double Rbsz = RZGP5*bsz;
        const double Rcsx = RZGP5*csx;
        const double Rcsy = RZGP5*csy;
        const double Rcsz = RZGP5*csz;
        const double Rdsx = RZGP5*dsx;
        const double Rdsy = RZGP5*dsy;
        const double Rdsz = RZGP5*dsz;
        const double Resx = RZGP5*esx;
        const double Resy = RZGP5*esy;
        const double Resz = RZGP5*esz;

        ERP[i][ 0] = Gbx*asa + Rysa + Rasy;
        ERP[i][ 1] = Gbx*asb + Rysb + Rasz;
        ERP[i][ 2] = Gbx*asc + Rysc;
        ERP[i][ 3] = Gbx*asd + Rysd + Rasx;
        ERP[i][ 4] = Gbx*ase + Ryse - SQR3I *Rasx;
        ERP[i][ 5] = Gby*asa + Rxsa + Rasx;
        ERP[i][ 6] = Gby*asb + Rxsb;
        ERP[i][ 7] = Gby*asc + Rxsc + Rasz;
        ERP[i][ 8] = Gby*asd + Rxsd - Rasy;
        ERP[i][ 9] = Gby*ase + Rxse - SQR3I *Rasy;
        ERP[i][10] = Gbz*asa;
        ERP[i][11] = Gbz*asb        + Rasx;
        ERP[i][12] = Gbz*asc        + Rasy;
        ERP[i][13] = Gbz*asd;
        ERP[i][14] = Gbz*ase        + SQR3I2*Rasz;

        ERP[i][15] = Gbx*bsa + Rzsa + Rbsy;
        ERP[i][16] = Gbx*bsb + Rzsb + Rbsz;
        ERP[i][17] = Gbx*bsc + Rzsc;
        ERP[i][18] = Gbx*bsd + Rzsd + Rbsx;
        ERP[i][19] = Gbx*bse + Rzse - SQR3I *Rbsx;
        ERP[i][20] = Gby*bsa        + Rbsx;
        ERP[i][21] = Gby*bsb;
        ERP[i][22] = Gby*bsc        + Rbsz;
        ERP[i][23] = Gby*bsd        - Rbsy;
        ERP[i][24] = Gby*bse        - SQR3I *Rbsy;
        ERP[i][25] = Gbz*bsa + Rxsa;
        ERP[i][26] = Gbz*bsb + Rxsb + Rbsx;
        ERP[i][27] = Gbz*bsc + Rxsc + Rbsy;
        ERP[i][28] = Gbz*bsd + Rxsd;
        ERP[i][29] = Gbz*bse + Rxse + SQR3I2*Rbsz;

        ERP[i][30] = Gbx*csa        + Rcsy;
        ERP[i][31] = Gbx*csb        + Rcsz;
        ERP[i][32] = Gbx*csc;
        ERP[i][33] = Gbx*csd        + Rcsx;
        ERP[i][34] = Gbx*cse        - SQR3I *Rcsx;
        ERP[i][35] = Gby*csa + Rzsa + Rcsx;
        ERP[i][36] = Gby*csb + Rzsb;
        ERP[i][37] = Gby*csc + Rzsc + Rcsz;
        ERP[i][38] = Gby*csd + Rzsd - Rcsy;
        ERP[i][39] = Gby*cse + Rzse - SQR3I *Rcsy;
        ERP[i][40] = Gbz*csa + Rysa;
        ERP[i][41] = Gbz*csb + Rysb + Rcsx;
        ERP[i][42] = Gbz*csc + Rysc + Rcsy;
        ERP[i][43] = Gbz*csd + Rysd;
        ERP[i][44] = Gbz*cse + Ryse + SQR3I2*Rcsz;

        ERP[i][45] = Gbx*dsa + Rxsa + Rdsy;
        ERP[i][46] = Gbx*dsb + Rxsb + Rdsz;
        ERP[i][47] = Gbx*dsc + Rxsc;
        ERP[i][48] = Gbx*dsd + Rxsd + Rdsx;
        ERP[i][49] = Gbx*dse + Rxse - SQR3I *Rdsx;
        ERP[i][50] = Gby*dsa - Rysa + Rdsx;
        ERP[i][51] = Gby*dsb - Rysb;
        ERP[i][52] = Gby*dsc - Rysc + Rdsz;
        ERP[i][53] = Gby*dsd - Rysd - Rdsy;
        ERP[i][54] = Gby*dse - Ryse - SQR3I *Rdsy;
        ERP[i][55] = Gbz*dsa;
        ERP[i][56] = Gbz*dsb        + Rdsx;
        ERP[i][57] = Gbz*dsc        + Rdsy;
        ERP[i][58] = Gbz*dsd;
        ERP[i][59] = Gbz*dse        + SQR3I2*Rdsz;

        ERP[i][60] = Gbx*esa - SQR3I *Rxsa + Resy;
        ERP[i][61] = Gbx*esb - SQR3I *Rxsb + Resz;
        ERP[i][62] = Gbx*esc - SQR3I *Rxsc;
        ERP[i][63] = Gbx*esd - SQR3I *Rxsd + Resx;
        ERP[i][64] = Gbx*ese - SQR3I *(Rxse + Resx);
        ERP[i][65] = Gby*esa - SQR3I *Rysa + Resx;
        ERP[i][66] = Gby*esb - SQR3I *Rysb;
        ERP[i][67] = Gby*esc - SQR3I *Rysc + Resz;
        ERP[i][68] = Gby*esd - SQR3I *Rysd - Resy;
        ERP[i][69] = Gby*ese - SQR3I *(Ryse + Resy);
        ERP[i][70] = Gbz*esa + SQR3I2*Rzsa;
        ERP[i][71] = Gbz*esb + SQR3I2*Rzsb + Resx;
        ERP[i][72] = Gbz*esc + SQR3I2*Rzsc + Resy;
        ERP[i][73] = Gbz*esd + SQR3I2*Rzsd;
        ERP[i][74] = Gbz*ese + SQR3I2*(Rzse + Resz);
    }
}


/*****************************************************************************/
/*                                       */
/*  Primitive 3-center OVP evaluation [ ABC ] for type [ DDD ]           */
/*                                       */
/*  np      ; number of first primitive pair data            */
/*  Zp      ; (Alpha+Beta) from orbital exponents Alpha, Beta    */
/*  px, py, pz  ; (Alpha*ax * Beta*bx) / (Alpha+Beta) , ...          */
/*  pax, pay, paz   ; ax, ay, az                         */
/*  pbx, pby, pbz   ; bx, by, bz                         */
/*  HP      ; Zeta**(-3/2) * exp(-(Alpha*Beta/Zeta) * (A-B)**2   */
/*  Gamma       ; orbital exponent Gamma for last GTO            */
/*  cx, cy, cz  ; coordinates of last GTO C              */
/*  cc      ; coefficient of last GTO C              */
/*  ERP     ; output primitive OVP                   */
/*          ERP[0,1,2,...,np-1][0] for [ DxyDxyDxy ]         */
/*          ERP[0,1,2,...,np-1][1] for [ DxyDxyDxz ]         */
/*          ERP[0,1,2,...,np-1][2] for [ DxyDxyDyz ]         */
/*          ERP[0,1,2,...,np-1][3] for [ DxyDxyDxx-yy ]      */
/*          ERP[0,1,2,...,np-1][4] for [ DxyDxyD3zz-rr ]         */
/*          ERP[0,1,2,...,np-1][5] for [ DxyDxzDxy ]         */
/*          ERP[0,1,2,...,np-1][6] for [ DxyDxzDxz ]         */
/*          ERP[0,1,2,...,np-1][7] for [ DxyDxzDyz ]         */
/*          ERP[0,1,2,...,np-1][8] for [ DxyDxzDxx-yy ]      */
/*          ERP[0,1,2,...,np-1][9] for [ DxyDxzD3zz-rr ]         */
/*                  .............                */
/*                                       */
/*****************************************************************************/
void DfOverlap::ovpDDD(const int np,
                       const double* px, const double* py, const double* pz,
                       const double* pax, const double* pay, const double* paz,
                       const double* pbx, const double* pby, const double* pbz,
                       const double Gamma,
                       const double cx, const double cy, const double cz,
                       const double cc,
                       const double* Zp, const double* HP, double** ERP)
{
    const double hc = M_PI*sqrt(M_PI)*cc;
    for (int i = 0; i < np; ++i) {
        const double RZG = 1.0 / (Zp[i] + Gamma);
        const double Roo = Zp[i]*Gamma*RZG;
        const double pcx = px[i] - cx;
        const double pcy = py[i] - cy;
        const double pcz = pz[i] - cz;
        const double t   = Roo * (pcx*pcx + pcy*pcy + pcz*pcz);
        const double sss = hc * HP[i] * RZG * sqrt(RZG) * exp(-t);

        const double Gx  = (Zp[i]*px[i] + Gamma*cx) * RZG;
        const double Gy  = (Zp[i]*py[i] + Gamma*cy) * RZG;
        const double Gz  = (Zp[i]*pz[i] + Gamma*cz) * RZG;
        const double Gax = Gx - pax[i];
        const double Gay = Gy - pay[i];
        const double Gaz = Gz - paz[i];
        const double Gbx = Gx - pbx[i];
        const double Gby = Gy - pby[i];
        const double Gbz = Gz - pbz[i];
        const double Gcx = Gx - cx;
        const double Gcy = Gy - cy;
        const double Gcz = Gz - cz;
        const double ABx = pax[i] - pbx[i];
        const double ABy = pay[i] - pby[i];
        const double ABz = paz[i] - pbz[i];
        // [ PSS ]
        const double RZGP5 = 0.5*RZG;
        const double xss   = Gax*sss;
        const double yss   = Gay*sss;
        const double zss   = Gaz*sss;
        // [ SPS ]
        const double sxs = xss + ABx*sss;
        const double sys = yss + ABy*sss;
        const double szs = zss + ABz*sss;
        // [ PPS ]
        const double Rsss = RZGP5*sss;
        const double xxs  = Gbx*xss + Rsss;
        const double xys  = Gby*xss;
        const double xzs  = Gbz*xss;
        const double yxs  = Gbx*yss;
        const double yys  = Gby*yss + Rsss;
        const double yzs  = Gbz*yss;
        const double zxs  = Gbx*zss;
        const double zys  = Gby*zss;
        const double zzs  = Gbz*zss + Rsss;
        // [ DSS ]  xy, xz, yz, xx-yy, 3zz-rr
        const double xxssww = Gax*xss;
        const double yyssww = Gay*yss;
        const double ass = Gax*yss;
        const double bss = Gax*zss;
        const double css = Gay*zss;
        const double dss = 0.5*(xxssww - yyssww);
        const double ess = SQR3I*(Gaz*zss
                                  - 0.5*(xxssww + yyssww));
        // [ DSS ] -2
        const double x2ss = xxssww  + Rsss;
        const double y2ss = yyssww  + Rsss;
        const double z2ss = Gaz*zss + Rsss;
        // [ SDS ]
        const double sas  = xys + ABx*sys;
        const double sbs  = xzs + ABx*szs;
        const double scs  = yzs + ABy*szs;
        const double sx2s = xxs + ABx*sxs;
        const double sy2s = yys + ABy*sys;
        const double sz2s = zzs + ABz*szs;
        const double sds  = 0.5*(sx2s - sy2s);
        const double ses  = SQR3I*(sz2s - 0.5*(sx2s + sy2s));
        // [ DPS ]
        const double Rxss = RZGP5*xss;
        const double Ryss = RZGP5*yss;
        const double Rzss = RZGP5*zss;
        const double axs = Gbx*ass + Ryss;
        const double ays = Gby*ass + Rxss;
        const double azs = Gbz*ass;
        const double bxs = Gbx*bss + Rzss;
        const double bys = Gby*bss;
        const double bzs = Gbz*bss + Rxss;
        const double cxs = Gbx*css;
        const double cys = Gby*css + Rzss;
        const double czs = Gbz*css + Ryss;
        const double dxs = Gbx*dss + Rxss;
        const double dys = Gby*dss - Ryss;
        const double dzs = Gbz*dss;
        const double exs = Gbx*ess - SQR3I *Rxss;
        const double eys = Gby*ess - SQR3I *Ryss;
        const double ezs = Gbz*ess + SQR3I2*Rzss;
        // [ DPS ] -2 (Dxx)xs
        const double x2xs = Gbx*x2ss + Rxss + Rxss;
        const double y2ys = Gby*y2ss + Ryss + Ryss;
        const double z2zs = Gbz*z2ss + Rzss + Rzss;
        // [ PDS ]
        const double xas = axs + ABy*xxs;
        const double xbs = bxs + ABz*xxs;
        const double xcs = azs + ABy*xzs;
        const double xds = 0.5*(x2xs + ABx*xxs - ays - ABy*xys);
        const double xes = SQR3I*(bzs + ABz*xzs
                                  -0.5*(x2xs + ABx*xxs + ays + ABy*xys));
        const double yas = ays + ABx*yys;
        const double ybs = cxs + ABz*yxs;
        const double ycs = cys + ABz*yys;
        const double yds = 0.5*(axs + ABx*yxs - y2ys - ABy*yys);
        const double yes = SQR3I*(czs + ABz*yzs
                                  -0.5*(axs + ABx*yxs + y2ys + ABy*yys));
        const double zas = bys + ABx*zys;
        const double zbs = bzs + ABx*zzs;
        const double zcs = czs + ABy*zzs;
        const double zds = 0.5*(bxs + ABx*zxs - cys - ABy*zys);
        const double zes = SQR3I*(z2zs + ABz*zzs
                                  -0.5*(bxs + ABx*zxs + cys + ABy*zys));
        // [ DPP ]
        const double Rass = RZGP5*ass;
        const double Rbss = RZGP5*bss;
        const double Rcss = RZGP5*css;
        const double Rdss = RZGP5*dss;
        const double Ress = RZGP5*ess;
        const double Rxxs = RZGP5*xxs;
        const double Rxys = RZGP5*xys;
        const double Rxzs = RZGP5*xzs;
        const double Ryxs = RZGP5*yxs;
        const double Ryys = RZGP5*yys;
        const double Ryzs = RZGP5*yzs;
        const double Rzxs = RZGP5*zxs;
        const double Rzys = RZGP5*zys;
        const double Rzzs = RZGP5*zzs;
        // a = xy
        const double axx = Gcx*axs + Ryxs + Rass;
        const double axy = Gcy*axs + Rxxs;
        const double axz = Gcz*axs;
        const double ayx = Gcx*ays + Ryys;
        const double ayy = Gcy*ays + Rxys + Rass;
        const double ayz = Gcz*ays;
        const double azx = Gcx*azs + Ryzs;
        const double azy = Gcy*azs + Rxzs;
        const double azz = Gcz*azs        + Rass;
        // b = xz
        const double bxx = Gcx*bxs + Rzxs + Rbss;
        const double bxy = Gcy*bxs;
        const double bxz = Gcz*bxs + Rxxs;
        const double byx = Gcx*bys + Rzys;
        const double byy = Gcy*bys        + Rbss;
        const double byz = Gcz*bys + Rxys;
        const double bzx = Gcx*bzs + Rzzs;
        const double bzy = Gcy*bzs;
        const double bzz = Gcz*bzs + Rxzs + Rbss;
        // c = yz
        const double cxx = Gcx*cxs        + Rcss;
        const double cxy = Gcy*cxs + Rzxs;
        const double cxz = Gcz*cxs + Ryxs;
        const double cyx = Gcx*cys;
        const double cyy = Gcy*cys + Rzys + Rcss;
        const double cyz = Gcz*cys + Ryys;
        const double czx = Gcx*czs;
        const double czy = Gcy*czs + Rzzs;
        const double czz = Gcz*czs + Ryzs + Rcss;
        // d = xx-yy
        const double dxx = Gcx*dxs + Rxxs + Rdss;
        const double dxy = Gcy*dxs - Ryxs;
        const double dxz = Gcz*dxs;
        const double dyx = Gcx*dys + Rxys;
        const double dyy = Gcy*dys - Ryys + Rdss;
        const double dyz = Gcz*dys;
        const double dzx = Gcx*dzs + Rxzs;
        const double dzy = Gcy*dzs - Ryzs;
        const double dzz = Gcz*dzs        + Rdss;
        // e = 3zz-rr
        const double exx = Gcx*exs - SQR3I *Rxxs + Ress;
        const double exy = Gcy*exs - SQR3I *Ryxs;
        const double exz = Gcz*exs + SQR3I2*Rzxs;
        const double eyx = Gcx*eys - SQR3I *Rxys;
        const double eyy = Gcy*eys - SQR3I *Ryys + Ress;
        const double eyz = Gcz*eys + SQR3I2*Rzys;
        const double ezx = Gcx*ezs - SQR3I *Rxzs;
        const double ezy = Gcy*ezs - SQR3I *Ryzs;
        const double ezz = Gcz*ezs + SQR3I2*Rzzs + Ress;
        // [ PDP ]
        const double Rsas = RZGP5*sas;
        const double Rsbs = RZGP5*sbs;
        const double Rscs = RZGP5*scs;
        const double Rsds = RZGP5*sds;
        const double Rses = RZGP5*ses;
        //const double Rxxs = RZGP5*xxs;
        //const double Rxys = RZGP5*xys;
        //const double Rxzs = RZGP5*xzs;
        //const double Ryxs = RZGP5*yxs;
        //const double Ryys = RZGP5*yys;
        //const double Ryzs = RZGP5*yzs;
        //const double Rzxs = RZGP5*zxs;
        //const double Rzys = RZGP5*zys;
        //const double Rzzs = RZGP5*zzs;
        const double xax  = Gcx*xas + Rsas + Rxys;
        const double xay  = Gcy*xas        + Rxxs;
        const double xaz  = Gcz*xas;
        const double xbx  = Gcx*xbs + Rsbs + Rxzs;
        const double xby  = Gcy*xbs;
        const double xbz  = Gcz*xbs        + Rxxs;
        const double xcx  = Gcx*xcs + Rscs;
        const double xcy  = Gcy*xcs        + Rxzs;
        const double xcz  = Gcz*xcs        + Rxys;
        const double xdx  = Gcx*xds + Rsds + Rxxs;
        const double xdy  = Gcy*xds        - Rxys;
        const double xdz  = Gcz*xds;
        const double xex  = Gcx*xes + Rses - SQR3I *Rxxs;
        const double xey  = Gcy*xes        - SQR3I *Rxys;
        const double xez  = Gcz*xes        + SQR3I2*Rxzs;
        const double yax  = Gcx*yas        + Ryys;
        const double yay  = Gcy*yas + Rsas + Ryxs;
        const double yaz  = Gcz*yas;
        const double ybx  = Gcx*ybs        + Ryzs;
        const double yby  = Gcy*ybs + Rsbs;
        const double ybz  = Gcz*ybs        + Ryxs;
        const double ycx  = Gcx*ycs;
        const double ycy  = Gcy*ycs + Rscs + Ryzs;
        const double ycz  = Gcz*ycs        + Ryys;
        const double ydx  = Gcx*yds        + Ryxs;
        const double ydy  = Gcy*yds + Rsds - Ryys;
        const double ydz  = Gcz*yds;
        const double yex  = Gcx*yes        - SQR3I *Ryxs;
        const double yey  = Gcy*yes + Rses - SQR3I *Ryys;
        const double yez  = Gcz*yes        + SQR3I2*Ryzs;
        const double zax  = Gcx*zas        + Rzys;
        const double zay  = Gcy*zas        + Rzxs;
        const double zaz  = Gcz*zas + Rsas;
        const double zbx  = Gcx*zbs        + Rzzs;
        const double zby  = Gcy*zbs;
        const double zbz  = Gcz*zbs + Rsbs + Rzxs;
        const double zcx  = Gcx*zcs;
        const double zcy  = Gcy*zcs        + Rzzs;
        const double zcz  = Gcz*zcs + Rscs + Rzys;
        const double zdx  = Gcx*zds        + Rzxs;
        const double zdy  = Gcy*zds        - Rzys;
        const double zdz  = Gcz*zds + Rsds;
        const double zex  = Gcx*zes        - SQR3I *Rzxs;
        const double zey  = Gcy*zes        - SQR3I *Rzys;
        const double zez  = Gcz*zes + Rses + SQR3I2*Rzzs;
        // [ DDS ]
        //const double Rxxs  = RZGP5*xxs;
        //const double Rxys  = RZGP5*xys;
        //const double Rxzs  = RZGP5*xzs;
        //const double Ryxs  = RZGP5*yxs;
        //const double Ryys  = RZGP5*yys;
        //const double Ryzs  = RZGP5*yzs;
        //const double Rzxs  = RZGP5*zxs;
        //const double Rzys  = RZGP5*zys;
        //const double Rzzs  = RZGP5*zzs;
        const double axxsw = Gbx*axs + Ryxs;
        const double ayysw = Gby*ays + Rxys;
        const double azzsw = Gbz*azs;
        const double bxxsw = Gbx*bxs + Rzxs;
        const double byysw = Gby*bys;
        const double bzzsw = Gbz*bzs + Rxzs;
        const double cxxsw = Gbx*cxs;
        const double cyysw = Gby*cys + Rzys;
        const double czzsw = Gbz*czs + Ryzs;
        const double dxxsw = Gbx*dxs + Rxxs;
        const double dyysw = Gby*dys - Ryys;
        const double dzzsw = Gbz*dzs;
        const double exxsw = Gbx*exs - SQR3I *Rxxs;
        const double eyysw = Gby*eys - SQR3I *Ryys;
        const double ezzsw = Gbz*ezs + SQR3I2*Rzzs;
        const double aas = Gbx*ays + Ryys;
        const double abs = Gbx*azs + Ryzs;
        const double acs = Gby*azs + Rxzs;
        const double ads = 0.5*(axxsw - ayysw);
        const double aes = SQR3I*(azzsw - 0.5*(axxsw + ayysw));
        const double bas = Gbx*bys + Rzys;
        const double bbs = Gbx*bzs + Rzzs;
        const double bcs = Gby*bzs;
        const double bds = 0.5*(bxxsw - byysw);
        const double bes = SQR3I*(bzzsw - 0.5*(bxxsw + byysw));
        const double cas = Gbx*cys;
        const double cbs = Gbx*czs;
        const double ccs = Gby*czs + Rzzs;
        const double cds = 0.5*(cxxsw - cyysw);
        const double ces = SQR3I*(czzsw - 0.5*(cxxsw + cyysw));
        const double das = Gbx*dys + Rxys;
        const double dbs = Gbx*dzs + Rxzs;
        const double dcs = Gby*dzs - Ryzs;
        const double dds = 0.5*(dxxsw - dyysw);
        const double des = SQR3I*(dzzsw - 0.5*(dxxsw + dyysw));
        const double eas = Gbx*eys - SQR3I*Rxys;
        const double ebs = Gbx*ezs - SQR3I*Rxzs;
        const double ecs = Gby*ezs - SQR3I*Ryzs;
        const double eds = 0.5*(exxsw - eyysw);
        const double ees = SQR3I*(ezzsw - 0.5*(exxsw + eyysw));
        // [ DDP ]
        const double Raxs = RZGP5*axs;
        const double Rays = RZGP5*ays;
        const double Razs = RZGP5*azs;
        const double Rbxs = RZGP5*bxs;
        const double Rbys = RZGP5*bys;
        const double Rbzs = RZGP5*bzs;
        const double Rcxs = RZGP5*cxs;
        const double Rcys = RZGP5*cys;
        const double Rczs = RZGP5*czs;
        const double Rdxs = RZGP5*dxs;
        const double Rdys = RZGP5*dys;
        const double Rdzs = RZGP5*dzs;
        const double Rexs = RZGP5*exs;
        const double Reys = RZGP5*eys;
        const double Rezs = RZGP5*ezs;
        const double Rxas = RZGP5*xas;
        const double Rxbs = RZGP5*xbs;
        const double Rxcs = RZGP5*xcs;
        const double Rxds = RZGP5*xds;
        const double Rxes = RZGP5*xes;
        const double Ryas = RZGP5*yas;
        const double Rybs = RZGP5*ybs;
        const double Rycs = RZGP5*ycs;
        const double Ryds = RZGP5*yds;
        const double Ryes = RZGP5*yes;
        const double Rzas = RZGP5*zas;
        const double Rzbs = RZGP5*zbs;
        const double Rzcs = RZGP5*zcs;
        const double Rzds = RZGP5*zds;
        const double Rzes = RZGP5*zes;

        const double aax = Gcx*aas + Ryas + Rays;
        const double aay = Gcy*aas + Rxas + Raxs;
        const double aaz = Gcz*aas;
        const double abx = Gcx*abs + Rybs + Razs;
        const double aby = Gcy*abs + Rxbs;
        const double abz = Gcz*abs        + Raxs;
        const double acx = Gcx*acs + Rycs;
        const double acy = Gcy*acs + Rxcs + Razs;
        const double acz = Gcz*acs        + Rays;
        const double adx = Gcx*ads + Ryds + Raxs;
        const double ady = Gcy*ads + Rxds - Rays;
        const double adz = Gcz*ads;
        const double aex = Gcx*aes + Ryes - SQR3I *Raxs;
        const double aey = Gcy*aes + Rxes - SQR3I *Rays;
        const double aez = Gcz*aes        + SQR3I2*Razs;

        const double bax = Gcx*bas + Rzas + Rbys;
        const double bay = Gcy*bas        + Rbxs;
        const double baz = Gcz*bas + Rxas;
        const double bbx = Gcx*bbs + Rzbs + Rbzs;
        const double bby = Gcy*bbs;
        const double bbz = Gcz*bbs + Rxbs + Rbxs;
        const double bcx = Gcx*bcs + Rzcs;
        const double bcy = Gcy*bcs        + Rbzs;
        const double bcz = Gcz*bcs + Rxcs + Rbys;
        const double bdx = Gcx*bds + Rzds + Rbxs;
        const double bdy = Gcy*bds        - Rbys;
        const double bdz = Gcz*bds + Rxds;
        const double bex = Gcx*bes + Rzes - SQR3I *Rbxs;
        const double bey = Gcy*bes        - SQR3I *Rbys;
        const double bez = Gcz*bes + Rxes + SQR3I2*Rbzs;

        const double cax = Gcx*cas        + Rcys;
        const double cay = Gcy*cas + Rzas + Rcxs;
        const double caz = Gcz*cas + Ryas;
        const double cbx = Gcx*cbs        + Rczs;
        const double cby = Gcy*cbs + Rzbs;
        const double cbz = Gcz*cbs + Rybs + Rcxs;
        const double ccx = Gcx*ccs;
        const double ccy = Gcy*ccs + Rzcs + Rczs;
        const double ccz = Gcz*ccs + Rycs + Rcys;
        const double cdx = Gcx*cds        + Rcxs;
        const double cdy = Gcy*cds + Rzds - Rcys;
        const double cdz = Gcz*cds + Ryds;
        const double cex = Gcx*ces        - SQR3I *Rcxs;
        const double cey = Gcy*ces + Rzes - SQR3I *Rcys;
        const double cez = Gcz*ces + Ryes + SQR3I2*Rczs;

        const double dax = Gcx*das + Rxas + Rdys;
        const double day = Gcy*das - Ryas + Rdxs;
        const double daz = Gcz*das;
        const double dbx = Gcx*dbs + Rxbs + Rdzs;
        const double dby = Gcy*dbs - Rybs;
        const double dbz = Gcz*dbs        + Rdxs;
        const double dcx = Gcx*dcs + Rxcs;
        const double dcy = Gcy*dcs - Rycs + Rdzs;
        const double dcz = Gcz*dcs        + Rdys;
        const double ddx = Gcx*dds + Rxds + Rdxs;
        const double ddy = Gcy*dds - Ryds - Rdys;
        const double ddz = Gcz*dds;
        const double dex = Gcx*des + Rxes - SQR3I *Rdxs;
        const double dey = Gcy*des - Ryes - SQR3I *Rdys;
        const double dez = Gcz*des        + SQR3I2*Rdzs;

        const double eax = Gcx*eas - SQR3I *Rxas + Reys;
        const double eay = Gcy*eas - SQR3I *Ryas + Rexs;
        const double eaz = Gcz*eas + SQR3I2*Rzas;
        const double ebx = Gcx*ebs - SQR3I *Rxbs + Rezs;
        const double eby = Gcy*ebs - SQR3I *Rybs;
        const double ebz = Gcz*ebs + SQR3I2*Rzbs + Rexs;
        const double ecx = Gcx*ecs - SQR3I *Rxcs;
        const double ecy = Gcy*ecs - SQR3I *Rycs + Rezs;
        const double ecz = Gcz*ecs + SQR3I2*Rzcs + Reys;
        const double edx = Gcx*eds - SQR3I *Rxds + Rexs;
        const double edy = Gcy*eds - SQR3I *Ryds - Reys;
        const double edz = Gcz*eds + SQR3I2*Rzds;
        const double eex = Gcx*ees - SQR3I *(Rxes + Rexs);
        const double eey = Gcy*ees - SQR3I *(Ryes + Reys);
        const double eez = Gcz*ees + SQR3I2*(Rzes + Rezs);
        // [ DDD ]
        const double Rxax = RZGP5*xax;
        const double Rxay = RZGP5*xay;
        const double Rxaz = RZGP5*xaz;
        const double Rxbx = RZGP5*xbx;
        const double Rxby = RZGP5*xby;
        const double Rxbz = RZGP5*xbz;
        const double Rxcx = RZGP5*xcx;
        const double Rxcy = RZGP5*xcy;
        const double Rxcz = RZGP5*xcz;
        const double Rxdx = RZGP5*xdx;
        const double Rxdy = RZGP5*xdy;
        const double Rxdz = RZGP5*xdz;
        const double Rxex = RZGP5*xex;
        const double Rxey = RZGP5*xey;
        const double Rxez = RZGP5*xez;
        const double Ryax = RZGP5*yax;
        const double Ryay = RZGP5*yay;
        const double Ryaz = RZGP5*yaz;
        const double Rybx = RZGP5*ybx;
        const double Ryby = RZGP5*yby;
        const double Rybz = RZGP5*ybz;
        const double Rycx = RZGP5*ycx;
        const double Rycy = RZGP5*ycy;
        const double Rycz = RZGP5*ycz;
        const double Rydx = RZGP5*ydx;
        const double Rydy = RZGP5*ydy;
        const double Rydz = RZGP5*ydz;
        const double Ryex = RZGP5*yex;
        const double Ryey = RZGP5*yey;
        const double Ryez = RZGP5*yez;
        const double Rzax = RZGP5*zax;
        const double Rzay = RZGP5*zay;
        const double Rzaz = RZGP5*zaz;
        const double Rzbx = RZGP5*zbx;
        const double Rzby = RZGP5*zby;
        const double Rzbz = RZGP5*zbz;
        const double Rzcx = RZGP5*zcx;
        const double Rzcy = RZGP5*zcy;
        const double Rzcz = RZGP5*zcz;
        const double Rzdx = RZGP5*zdx;
        const double Rzdy = RZGP5*zdy;
        const double Rzdz = RZGP5*zdz;
        const double Rzex = RZGP5*zex;
        const double Rzey = RZGP5*zey;
        const double Rzez = RZGP5*zez;

        const double Raxx = RZGP5*axx;
        const double Raxy = RZGP5*axy;
        const double Raxz = RZGP5*axz;
        const double Rayx = RZGP5*ayx;
        const double Rayy = RZGP5*ayy;
        const double Rayz = RZGP5*ayz;
        const double Razx = RZGP5*azx;
        const double Razy = RZGP5*azy;
        const double Razz = RZGP5*azz;
        const double Rbxx = RZGP5*bxx;
        const double Rbxy = RZGP5*bxy;
        const double Rbxz = RZGP5*bxz;
        const double Rbyx = RZGP5*byx;
        const double Rbyy = RZGP5*byy;
        const double Rbyz = RZGP5*byz;
        const double Rbzx = RZGP5*bzx;
        const double Rbzy = RZGP5*bzy;
        const double Rbzz = RZGP5*bzz;
        const double Rcxx = RZGP5*cxx;
        const double Rcxy = RZGP5*cxy;
        const double Rcxz = RZGP5*cxz;
        const double Rcyx = RZGP5*cyx;
        const double Rcyy = RZGP5*cyy;
        const double Rcyz = RZGP5*cyz;
        const double Rczx = RZGP5*czx;
        const double Rczy = RZGP5*czy;
        const double Rczz = RZGP5*czz;
        const double Rdxx = RZGP5*dxx;
        const double Rdxy = RZGP5*dxy;
        const double Rdxz = RZGP5*dxz;
        const double Rdyx = RZGP5*dyx;
        const double Rdyy = RZGP5*dyy;
        const double Rdyz = RZGP5*dyz;
        const double Rdzx = RZGP5*dzx;
        const double Rdzy = RZGP5*dzy;
        const double Rdzz = RZGP5*dzz;
        const double Rexx = RZGP5*exx;
        const double Rexy = RZGP5*exy;
        const double Rexz = RZGP5*exz;
        const double Reyx = RZGP5*eyx;
        const double Reyy = RZGP5*eyy;
        const double Reyz = RZGP5*eyz;
        const double Rezx = RZGP5*ezx;
        const double Rezy = RZGP5*ezy;
        const double Rezz = RZGP5*ezz;

        const double aaxxw = Gcx*aax + Ryax + Rayx;
        const double aayyw = Gcy*aay + Rxay + Raxy;
        const double aazzw = Gcz*aaz;
        const double abxxw = Gcx*abx + Rybx + Razx;
        const double abyyw = Gcy*aby + Rxby;
        const double abzzw = Gcz*abz         + Raxz;
        const double acxxw = Gcx*acx + Rycx;
        const double acyyw = Gcy*acy + Rxcy + Razy;
        const double aczzw = Gcz*acz         + Rayz;
        const double adxxw = Gcx*adx + Rydx + Raxx;
        const double adyyw = Gcy*ady + Rxdy - Rayy;
        const double adzzw = Gcz*adz;
        const double aexxw = Gcx*aex + Ryex - SQR3I *Raxx;
        const double aeyyw = Gcy*aey + Rxey - SQR3I *Rayy;
        const double aezzw = Gcz*aez         + SQR3I2*Razz;

        const double baxxw = Gcx*bax + Rzax + Rbyx;
        const double bayyw = Gcy*bay         + Rbxy;
        const double bazzw = Gcz*baz + Rxaz;
        const double bbxxw = Gcx*bbx + Rzbx + Rbzx;
        const double bbyyw = Gcy*bby;
        const double bbzzw = Gcz*bbz + Rxbz + Rbxz;
        const double bcxxw = Gcx*bcx + Rzcx;
        const double bcyyw = Gcy*bcy         + Rbzy;
        const double bczzw = Gcz*bcz + Rxcz + Rbyz;
        const double bdxxw = Gcx*bdx + Rzdx + Rbxx;
        const double bdyyw = Gcy*bdy         - Rbyy;
        const double bdzzw = Gcz*bdz + Rxdz;
        const double bexxw = Gcx*bex + Rzex - SQR3I *Rbxx;
        const double beyyw = Gcy*bey         - SQR3I *Rbyy;
        const double bezzw = Gcz*bez + Rxez + SQR3I2*Rbzz;

        const double caxxw = Gcx*cax         + Rcyx;
        const double cayyw = Gcy*cay + Rzay + Rcxy;
        const double cazzw = Gcz*caz + Ryaz;
        const double cbxxw = Gcx*cbx         + Rczx;
        const double cbyyw = Gcy*cby + Rzby;
        const double cbzzw = Gcz*cbz + Rybz + Rcxz;
        const double ccxxw = Gcx*ccx;
        const double ccyyw = Gcy*ccy + Rzcy + Rczy;
        const double cczzw = Gcz*ccz + Rycz + Rcyz;
        const double cdxxw = Gcx*cdx         + Rcxx;
        const double cdyyw = Gcy*cdy + Rzdy - Rcyy;
        const double cdzzw = Gcz*cdz + Rydz;
        const double cexxw = Gcx*cex         - SQR3I *Rcxx;
        const double ceyyw = Gcy*cey + Rzey - SQR3I *Rcyy;
        const double cezzw = Gcz*cez + Ryez + SQR3I2*Rczz;

        const double daxxw = Gcx*dax + Rxax + Rdyx;
        const double dayyw = Gcy*day - Ryay + Rdxy;
        const double dazzw = Gcz*daz;
        const double dbxxw = Gcx*dbx + Rxbx + Rdzx;
        const double dbyyw = Gcy*dby - Ryby;
        const double dbzzw = Gcz*dbz         + Rdxz;
        const double dcxxw = Gcx*dcx + Rxcx;
        const double dcyyw = Gcy*dcy - Rycy + Rdzy;
        const double dczzw = Gcz*dcz         + Rdyz;
        const double ddxxw = Gcx*ddx + Rxdx + Rdxx;
        const double ddyyw = Gcy*ddy - Rydy - Rdyy;
        const double ddzzw = Gcz*ddz;
        const double dexxw = Gcx*dex + Rxex - SQR3I *Rdxx;
        const double deyyw = Gcy*dey - Ryey - SQR3I *Rdyy;
        const double dezzw = Gcz*dez         + SQR3I2*Rdzz;

        const double eaxxw = Gcx*eax - SQR3I *Rxax + Reyx;
        const double eayyw = Gcy*eay - SQR3I *Ryay + Rexy;
        const double eazzw = Gcz*eaz + SQR3I2*Rzaz;
        const double ebxxw = Gcx*ebx - SQR3I *Rxbx + Rezx;
        const double ebyyw = Gcy*eby - SQR3I *Ryby;
        const double ebzzw = Gcz*ebz + SQR3I2*Rzbz + Rexz;
        const double ecxxw = Gcx*ecx - SQR3I *Rxcx;
        const double ecyyw = Gcy*ecy - SQR3I *Rycy + Rezy;
        const double eczzw = Gcz*ecz + SQR3I2*Rzcz + Reyz;
        const double edxxw = Gcx*edx - SQR3I *Rxdx + Rexx;
        const double edyyw = Gcy*edy - SQR3I *Rydy - Reyy;
        const double edzzw = Gcz*edz + SQR3I2*Rzdz;
        const double eexxw = Gcx*eex - SQR3I *(Rxex + Rexx);
        const double eeyyw = Gcy*eey - SQR3I *(Ryey + Reyy);
        const double eezzw = Gcz*eez + SQR3I2*(Rzez + Rezz);

        ERP[i][  0] = Gcx*aay + Ryay + Rayy;
        ERP[i][  1] = Gcz*aax;
        ERP[i][  2] = Gcz*aay;
        ERP[i][  3] = 0.5*(aaxxw - aayyw);
        ERP[i][  4] = SQR3I*(aazzw - 0.5*(aaxxw + aayyw));
        ERP[i][  5] = Gcx*aby + Ryby + Razy;
        ERP[i][  6] = Gcz*abx        + Raxx;
        ERP[i][  7] = Gcz*aby        + Raxy;
        ERP[i][  8] = 0.5*(abxxw - abyyw);
        ERP[i][  9] = SQR3I*(abzzw - 0.5*(abxxw + abyyw));
        ERP[i][ 10] = Gcx*acy + Rycy;
        ERP[i][ 11] = Gcz*acx        + Rayx;
        ERP[i][ 12] = Gcz*acy        + Rayy;
        ERP[i][ 13] = 0.5*(acxxw - acyyw);
        ERP[i][ 14] = SQR3I*(aczzw - 0.5*(acxxw + acyyw));
        ERP[i][ 15] = Gcx*ady + Rydy + Raxy;
        ERP[i][ 16] = Gcz*adx;
        ERP[i][ 17] = Gcz*ady;
        ERP[i][ 18] = 0.5*(adxxw - adyyw);
        ERP[i][ 19] = SQR3I*(adzzw - 0.5*(adxxw + adyyw));
        ERP[i][ 20] = Gcx*aey + Ryey - SQR3I *Raxy;
        ERP[i][ 21] = Gcz*aex        + SQR3I2*Razx;
        ERP[i][ 22] = Gcz*aey        + SQR3I2*Razy;
        ERP[i][ 23] = 0.5*(aexxw - aeyyw);
        ERP[i][ 24] = SQR3I*(aezzw - 0.5*(aexxw + aeyyw));

        ERP[i][ 25] = Gcx*bay + Rzay + Rbyy;
        ERP[i][ 26] = Gcz*bax + Rxax;
        ERP[i][ 27] = Gcz*bay + Rxay;
        ERP[i][ 28] = 0.5*(baxxw - bayyw);
        ERP[i][ 29] = SQR3I*(bazzw - 0.5*(baxxw + bayyw));
        ERP[i][ 30] = Gcx*bby + Rzby + Rbzy;
        ERP[i][ 31] = Gcz*bbx + Rxbx + Rbxx;
        ERP[i][ 32] = Gcz*bby + Rxby + Rbxy;
        ERP[i][ 33] = 0.5*(bbxxw - bbyyw);
        ERP[i][ 34] = SQR3I*(bbzzw - 0.5*(bbxxw + bbyyw));
        ERP[i][ 35] = Gcx*bcy + Rzcy;
        ERP[i][ 36] = Gcz*bcx + Rxcx + Rbyx;
        ERP[i][ 37] = Gcz*bcy + Rxcy + Rbyy;
        ERP[i][ 38] = 0.5*(bcxxw - bcyyw);
        ERP[i][ 39] = SQR3I*(bczzw - 0.5*(bcxxw + bcyyw));
        ERP[i][ 40] = Gcx*bdy + Rzdy + Rbxy;
        ERP[i][ 41] = Gcz*bdx + Rxdx;
        ERP[i][ 42] = Gcz*bdy + Rxdy;
        ERP[i][ 43] = 0.5*(bdxxw - bdyyw);
        ERP[i][ 44] = SQR3I*(bdzzw - 0.5*(bdxxw + bdyyw));
        ERP[i][ 45] = Gcx*bey + Rzey - SQR3I *Rbxy;
        ERP[i][ 46] = Gcz*bex + Rxex + SQR3I2*Rbzx;
        ERP[i][ 47] = Gcz*bey + Rxey + SQR3I2*Rbzy;
        ERP[i][ 48] = 0.5*(bexxw - beyyw);
        ERP[i][ 49] = SQR3I*(bezzw - 0.5*(bexxw + beyyw));

        ERP[i][ 50] = Gcx*cay        + Rcyy;
        ERP[i][ 51] = Gcz*cax + Ryax;
        ERP[i][ 52] = Gcz*cay + Ryay;
        ERP[i][ 53] = 0.5*(caxxw - cayyw);
        ERP[i][ 54] = SQR3I*(cazzw - 0.5*(caxxw + cayyw));
        ERP[i][ 55] = Gcx*cby        + Rczy;
        ERP[i][ 56] = Gcz*cbx + Rybx + Rcxx;
        ERP[i][ 57] = Gcz*cby + Ryby + Rcxy;
        ERP[i][ 58] = 0.5*(cbxxw - cbyyw);
        ERP[i][ 59] = SQR3I*(cbzzw - 0.5*(cbxxw + cbyyw));
        ERP[i][ 60] = Gcx*ccy;
        ERP[i][ 61] = Gcz*ccx + Rycx + Rcyx;
        ERP[i][ 62] = Gcz*ccy + Rycy + Rcyy;
        ERP[i][ 63] = 0.5*(ccxxw - ccyyw);
        ERP[i][ 64] = SQR3I*(cczzw - 0.5*(ccxxw + ccyyw));
        ERP[i][ 65] = Gcx*cdy        + Rcxy;
        ERP[i][ 66] = Gcz*cdx + Rydx;
        ERP[i][ 67] = Gcz*cdy + Rydy;
        ERP[i][ 68] = 0.5*(cdxxw - cdyyw);
        ERP[i][ 69] = SQR3I*(cdzzw - 0.5*(cdxxw + cdyyw));
        ERP[i][ 70] = Gcx*cey        - SQR3I *Rcxy;
        ERP[i][ 71] = Gcz*cex + Ryex + SQR3I2*Rczx;
        ERP[i][ 72] = Gcz*cey + Ryey + SQR3I2*Rczy;
        ERP[i][ 73] = 0.5*(cexxw - ceyyw);
        ERP[i][ 74] = SQR3I*(cezzw - 0.5*(cexxw + ceyyw));

        ERP[i][ 75] = Gcx*day + Rxay + Rdyy;
        ERP[i][ 76] = Gcz*dax;
        ERP[i][ 77] = Gcz*day;
        ERP[i][ 78] = 0.5*(daxxw - dayyw);
        ERP[i][ 79] = SQR3I*(dazzw - 0.5*(daxxw + dayyw));
        ERP[i][ 80] = Gcx*dby + Rxby + Rdzy;
        ERP[i][ 81] = Gcz*dbx        + Rdxx;
        ERP[i][ 82] = Gcz*dby        + Rdxy;
        ERP[i][ 83] = 0.5*(dbxxw - dbyyw);
        ERP[i][ 84] = SQR3I*(dbzzw - 0.5*(dbxxw + dbyyw));
        ERP[i][ 85] = Gcx*dcy + Rxcy;
        ERP[i][ 86] = Gcz*dcx        + Rdyx;
        ERP[i][ 87] = Gcz*dcy        + Rdyy;
        ERP[i][ 88] = 0.5*(dcxxw - dcyyw);
        ERP[i][ 89] = SQR3I*(dczzw - 0.5*(dcxxw + dcyyw));
        ERP[i][ 90] = Gcx*ddy + Rxdy + Rdxy;
        ERP[i][ 91] = Gcz*ddx;
        ERP[i][ 92] = Gcz*ddy;
        ERP[i][ 93] = 0.5*(ddxxw - ddyyw);
        ERP[i][ 94] = SQR3I*(ddzzw - 0.5*(ddxxw + ddyyw));
        ERP[i][ 95] = Gcx*dey + Rxey - SQR3I *Rdxy;
        ERP[i][ 96] = Gcz*dex        + SQR3I2*Rdzx;
        ERP[i][ 97] = Gcz*dey        + SQR3I2*Rdzy;
        ERP[i][ 98] = 0.5*(dexxw - deyyw);
        ERP[i][ 99] = SQR3I*(dezzw - 0.5*(dexxw + deyyw));

        ERP[i][100] = Gcx*eay - SQR3I *Rxay + Reyy;
        ERP[i][101] = Gcz*eax + SQR3I2*Rzax;
        ERP[i][102] = Gcz*eay + SQR3I2*Rzay;
        ERP[i][103] = 0.5*(eaxxw - eayyw);
        ERP[i][104] = SQR3I*(eazzw - 0.5*(eaxxw + eayyw));
        ERP[i][105] = Gcx*eby - SQR3I *Rxby + Rezy;
        ERP[i][106] = Gcz*ebx + SQR3I2*Rzbx + Rexx;
        ERP[i][107] = Gcz*eby + SQR3I2*Rzby + Rexy;
        ERP[i][108] = 0.5*(ebxxw - ebyyw);
        ERP[i][109] = SQR3I*(ebzzw - 0.5*(ebxxw + ebyyw));
        ERP[i][110] = Gcx*ecy - SQR3I *Rxcy;
        ERP[i][111] = Gcz*ecx + SQR3I2*Rzcx + Reyx;
        ERP[i][112] = Gcz*ecy + SQR3I2*Rzcy + Reyy;
        ERP[i][113] = 0.5*(ecxxw - ecyyw);
        ERP[i][114] = SQR3I*(eczzw - 0.5*(ecxxw + ecyyw));
        ERP[i][115] = Gcx*edy - SQR3I *Rxdy + Rexy;
        ERP[i][116] = Gcz*edx + SQR3I2*Rzdx;
        ERP[i][117] = Gcz*edy + SQR3I2*Rzdy;
        ERP[i][118] = 0.5*(edxxw - edyyw);
        ERP[i][119] = SQR3I*(edzzw - 0.5*(edxxw + edyyw));
        ERP[i][120] = Gcx*eey - SQR3I *(Rxey + Rexy);
        ERP[i][121] = Gcz*eex + SQR3I2*(Rzex + Rezx);
        ERP[i][122] = Gcz*eey + SQR3I2*(Rzey + Rezy);
        ERP[i][123] = 0.5*(eexxw - eeyyw);
        ERP[i][124] = SQR3I*(eezzw - 0.5*(eexxw + eeyyw));
    }
}

