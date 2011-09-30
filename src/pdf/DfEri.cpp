#include <stdlib.h>
#include <cmath>
#include <cassert>

#include "DfEri.h"
#include "TlSymmetricMatrix.h"
#include "TlVector.h"
#include "TlFmt.h"
#include "TlLogX.h"
#include "TlTime.h"

#define MAX_SHELL (2) // s=0, p=1, d=2
#define DEFAULT_BLOCK_SIZE 1000

const double DfEri::FPAI   = 2.0*pow(M_PI, 2.5);
const double DfEri::D33    = 1.0/0.03;
const double DfEri::SQR3I  = 1.0/sqrt(3.0);
const double DfEri::SQR3I2 = 2.0/sqrt(3.0);
const int DfEri::MAX_TYPE = 3; // s, p, d

DfEri::DfEri(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam), pOrbitalInfo_(NULL), m_dStoredCutValue(0.0)
{
    const TlSerializeData& pdfParam = *pPdfParam;
    this->m_dStoredCutValue = pdfParam["cut-value"].getDouble();
    this->m_dCutValue = m_dStoredCutValue;

    this->blockSize_ = DEFAULT_BLOCK_SIZE;
    if (pdfParam["block-size"].getStr().empty() != true) {
        this->blockSize_ = pdfParam["block-size"].getInt();
    }

    this->pOrbitalInfo_ = new TlOrbitalInfo((*pPdfParam)["coordinates"],
                                            (*pPdfParam)["basis_set"]);
    this->pOrbitalInfo_Density_ = new TlOrbitalInfo_Density((*pPdfParam)["coordinates"],
                                                            (*pPdfParam)["basis_set_auxD"]);
    
    this->getMemory();
    this->initialize();
}


DfEri::~DfEri()
{
    if (this->pOrbitalInfo_ != NULL) {
        delete this->pOrbitalInfo_;
        this->pOrbitalInfo_ = NULL;
    }
    if (this->pOrbitalInfo_Density_ != NULL) {
        delete this->pOrbitalInfo_Density_;
        this->pOrbitalInfo_Density_ = NULL;
    }

    for (int i = 0; i < 60; ++i) {
        delete[] *(ADAT+i);
    }
    delete[] ADAT;
    for (int i = 0; i < 7; ++i) {
        delete[] *(FDAT+i);
    }
    delete[] FDAT;
    delete[] SDAT;
    delete[] EDAT;
    delete[] GA;
    delete[] RMI;
    delete[] TF;
}


void DfEri::getMemory()
{
    TF = new double[25];
    for (int i = 0; i < 25; i++) {
        TF[i] = 0.0;
    }

    RMI = new double[24];
    for (int i = 0; i < 24; i++) {
        RMI[i] = 0.0;
    }

    GA = new double[6];
    for (int i = 0; i < 6; i++) {
        GA[i] = 0.0;
    }

    EDAT = new double[ 1901 ];
    for (int i = 0; i < 1901; i++) {
        EDAT[i] = 0.0;
    }

    SDAT = new double[ 1901 ];
    for (int i = 0; i < 1901; i++) {
        SDAT[i] = 0.0;
    }

    FDAT = new double*[ 7 ];
    for (int i = 0; i < 7; ++i) {
        *(FDAT+i) = new double[ 1901 ];
        for (int j = 0; j < 1901; j++) {
            FDAT[i][j] = 0.0;
        }
    }

    ADAT = new double*[ 60 ];
    for (int i = 0; i < 60; ++i) {
        *(ADAT+i) = new double[ 1901 ];
        for (int j = 0; j < 1901; j++) {
            ADAT[i][j] = 0.0;
        }
    }
}


void DfEri::initialize()
{
    this->makeTable();
    this->auxSet();
}


// 統計用
void DfEri::initializeCounter()
{
    this->TcountSS = 0;
    this->TcountPS = 0;
    this->TcountPP = 0;
    this->TcountDS = 0;
    this->TcountDP = 0;
    this->TcountDD = 0;
    this->CcountSS = 0;
    this->CcountPS = 0;
    this->CcountPP = 0;
    this->CcountDS = 0;
    this->CcountDP = 0;
    this->CcountDD = 0;
}


void DfEri::countupTotal(int ity, int jty)
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


void DfEri::countupCutoff(int ity, int jty)
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


void DfEri::cutoffReport()
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


void DfEri::cutoffReport(const std::string& shell, const std::size_t cutoffCount, const std::size_t totalCount)
{
    double rate = 0.0;
    if (totalCount > 0) {
        rate = (double)cutoffCount / (double)totalCount * 100.0;
    }
    std::stringstream ss;
    ss << " " << shell << " shell "
       << cutoffCount << " / "
       << totalCount << " "
       << TlUtils::format(" (%6.2f%%)\n", rate);
    this->logger(ss.str());
}


void DfEri::getSab(TlSymmetricMatrix* pSab)
{
    assert(pSab != NULL);
    pSab->resize(this->m_nNumOfAux);

    this->getSab_core(pSab);

    this->loggerTime(" finalize");
    // nothing to do
}


void DfEri::getSab_core(TlMatrixObject* pSab)
{
    assert(pSab != NULL);

    this->loggerTime(" integral");
    const std::size_t maxNumOfPairs = this->blockSize_;

    this->getQueue(NULL, this->shellList_Dens_,
                   maxNumOfPairs, true); // initialize only

    std::vector<IJShellPair> ij;
    do {
        ij = this->getQueue(NULL, this->shellList_Dens_, maxNumOfPairs);
        this->sabcalc(ij, pSab);
    } while (ij.empty() != true);
}


void DfEri::getDeltaT(const TlSymmetricMatrix& deltaPpq, TlVector* pDeltaT)
{
    // Set new cutvalue
    {
        const double MAXdeltaPpq = deltaPpq.getMaxAbsoluteElement();
        if ((MAXdeltaPpq > this->m_dCutValue) && (MAXdeltaPpq < 1.0)) {
            this->m_dCutValue /= MAXdeltaPpq;
            this->logger(TlUtils::format("new cutvalue in DfEri.getdeltaT is %.2e\n", this->m_dCutValue));
        } else {
            this->m_dCutValue = this->m_dStoredCutValue;
        }
    }

    this->getDeltaT_core(&deltaPpq, pDeltaT);
}

void DfEri::getDeltaT_core(const TlMatrixObject* deltaPpq, TlVectorObject* pDeltaT)
{
    this->loggerTime(" ERI start");

    TcountSS = TcountPS = TcountPP = TcountDS = TcountDP = TcountDD = 0;
    CcountSS = CcountPS = CcountPP = CcountDS = CcountDP = CcountDD = 0;

    std::list<IJShellPair> IJShellPairList;

    const std::size_t maxNumOfPairs = this->blockSize_;

    this->getQueue(this->pOrbitalInfo_, this->shellList_,
                   maxNumOfPairs, true); // initialize only
    std::vector<IJShellPair> ij; 
    do {
        ij = this->getQueue(this->pOrbitalInfo_, this->shellList_, maxNumOfPairs);
        this->ericalcDT(ij, deltaPpq, pDeltaT);
    } while (ij.empty() != true);

    this->cutoffReport();

    this->loggerTime(" ERI end");
}

// delta [pq|alpha]
void DfEri::getdeltaHpqA(const TlVector& deltaRho, TlSymmetricMatrix& deltaH)
{
    this->loggerTime(" 'deltaRho * <pq|alpha>' start");

    // Set new cutvalue
    const double MAXdeltaRho = deltaRho.getMaxAbsoluteElement();
    this->m_dCutValue = this->m_dStoredCutValue;
    if ((MAXdeltaRho > this->m_dCutValue) && (MAXdeltaRho < 1.0)) {
        this->m_dCutValue /= MAXdeltaRho;
        this->logger(TlUtils::format("new cutoff value of (deltaRho * <pq|alpha>) is %.2e\n", this->m_dCutValue));
    }

    TcountSS = TcountPS = TcountPP = TcountDS = TcountDP = TcountDD = 0;
    CcountSS = CcountPS = CcountPP = CcountDS = CcountDP = CcountDD = 0;

    const std::size_t maxNumOfPairs = this->blockSize_;

    this->getQueue(this->pOrbitalInfo_, this->shellList_,
                   maxNumOfPairs, true); // initialize only

    std::vector<IJShellPair> ij;
    do {
        ij = this->getQueue(this->pOrbitalInfo_, this->shellList_, maxNumOfPairs);
        this->ericalcDH(ij, deltaRho, &deltaH);
    } while (ij.empty() != true);

    this->finalizeIntegral(deltaH);
    this->cutoffReport();

    this->loggerTime(" end");
}


void DfEri::getDeltaHpqAForEri2(const TlVector& deltaRho, TlSymmetricMatrix& deltaH, const double dCutoffCoef)
{
    assert(0 < dCutoffCoef);

    // START
    this->loggerTime(" start");

    const std::size_t maxNumOfPairs = this->blockSize_;
    this->getQueue(this->pOrbitalInfo_, this->shellList_,
                   maxNumOfPairs, true); // initialize only

    std::vector<IJShellPair> ij;
    do {
        ij = this->getQueue(this->pOrbitalInfo_, this->shellList_,
                            maxNumOfPairs);
        this->ericalcDH(ij, deltaRho, &deltaH);
    } while (ij.empty() != true);

    this->finalizeIntegral(deltaH);

    // finish
    this->loggerTime(" end");
}


std::vector<DfObject::IJShellPair>
DfEri::getQueue(const TlOrbitalInfo* pOrbitalInfo,
                const std::vector<std::vector<std::size_t> >& shellList,
                const int maxNumOfPairs, const bool initialize)
{
    std::vector<IJShellPair> ijShellPairList;

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

        return ijShellPairList;
    }

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
                        if (this->isCutoffUsingSchwartzInequality(*pOrbitalInfo,
                                                                  index_i, index_j,
                                                                  this->m_dCutValue) == true) {
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


bool DfEri::getQueueX(const std::vector<std::vector<std::size_t> >& shellList,
                      int* pShellType, index_type* pShell,
                      const bool initialize)
{
    static int iShellType;
    static int jShellType;
    static index_type iShellIndex;
    
    if (initialize == true) {
        iShellType = MAX_TYPE -1;
        jShellType = iShellType;
        iShellIndex = 0;
        
        *pShellType = -1;
        *pShell = -1;
        
        return true;
    }

    while (iShellType >= 0) {
        while (jShellType >= 0) {
            assert(iShellType >= jShellType);
            *pShellType = iShellType * MAX_TYPE + jShellType;
            while ((std::size_t)iShellIndex < shellList[iShellType].size()) {
                *pShell = shellList[iShellType][iShellIndex];

                ++iShellIndex;
                return true;
            }
            --jShellType;
            iShellIndex = 0;
        }
        --iShellType;
        jShellType = iShellType;
    }

    return false;
}


void DfEri::finalizeIntegral(TlSymmetricMatrix& rMatrix)
{
    // do nothing
    //std::cerr << "DfEri::finalizeIntegral(m)" << std::endl;

}

void DfEri::finalizeIntegral(TlVector& rVector)
{
    // do nothing
    //std::cerr << "DfEri::finalizeIntegral(v)" << std::endl;
}

void DfEri::sabcalc(const std::vector<IJShellPair>& IJShellList, TlMatrixObject* Sab)
{
    const double dCutValue = this->m_dCutValue;
    const double one = 1.0;
    const double zero = 0.0;
    double O[3];
    O[0] = 0.0;
    O[1] = 0.0;
    O[2] = 0.0;

    const int maxNumOfPairs = IJShellList.size();
#pragma omp parallel
    {
        double ERI[(2 *MAX_SHELL+1) *(2 *MAX_SHELL+1) *(2 *MAX_SHELL+1)];
#pragma omp for
        for (int pairIndex = 0; pairIndex < maxNumOfPairs; ++pairIndex) {
            const std::size_t orbA = IJShellList[pairIndex].nIShell;
            const std::size_t orbB = IJShellList[pairIndex].nJShell;

            const int nqA = this->pOrbitalInfo_Density_->getShellType(orbA);
            const int iAngA = nqA * 2 + 1;
            const double coefA = this->pOrbitalInfo_Density_->getCoefficient(orbA, 0);
            const double zetaA = this->pOrbitalInfo_Density_->getExponent(orbA, 0);
            const TlPosition posA = this->pOrbitalInfo_Density_->getPosition(orbA);
            double A[3];
            A[0] = posA.x();
            A[1] = posA.y();
            A[2] = posA.z();

            const int nqB = this->pOrbitalInfo_Density_->getShellType(orbB);
            const int iAngB = nqB * 2 + 1;
            const double coefB = this->pOrbitalInfo_Density_->getCoefficient(orbB, 0);
            const double zetaB = this->pOrbitalInfo_Density_->getExponent(orbB, 0);
            const TlPosition posB = this->pOrbitalInfo_Density_->getPosition(orbB);
            double B[3];
            B[0] = posB.x();
            B[1] = posB.y();
            B[2] = posB.z();

            this->ericalc(nqA, 0, nqB,
                          1, &coefA, &zetaA, A,
                          1, &one, &zero, O,
                          coefB, zetaB, B, ERI);

#pragma omp critical (DfEri_sabcalc)
            {
                int iww = 0;
                for (int iAng = 0; iAng < iAngA; ++iAng) {
                    const std::size_t I = orbA + iAng;

                    for (int jAng = 0; jAng < iAngB; ++jAng) {
                        const std::size_t J = orbB + jAng;

                        if ((orbA != orbB) || (I >= J)) {
                            if (std::fabs(ERI[iww]) >= dCutValue) {
                                Sab->add(I, J, ERI[iww]);
                            }
                        }
                        ++iww;
                    }
                }
            }
        }
    }
}


void DfEri::ericalcDT(const std::vector<IJShellPair>& IJShellList,
                      const TlMatrixObject* deltaPpq, TlVectorObject* pDeltaT)
{
    const double dCutValue = this->m_dCutValue;
    const std::size_t numOfAuxDens = this->m_nNumOfAux;

    const int maxNumOfPairs = IJShellList.size();
#pragma omp parallel
    {
        double ERI[(2 *MAX_SHELL+1) *(2 *MAX_SHELL+1) *(2 *MAX_SHELL+1)];
#pragma omp for
        for (int pairIndex = 0; pairIndex < maxNumOfPairs; ++pairIndex) {
            const std::size_t orbA = IJShellList[pairIndex].nIShell;
            const std::size_t orbB = IJShellList[pairIndex].nJShell;

            const int nqA = this->pOrbitalInfo_->getShellType(orbA);
            const int iAngA  = nqA * 2 + 1;
            const int numOfContractionA = this->pOrbitalInfo_->getCgtoContraction(orbA);
            double* pCoefA = new double[numOfContractionA];
            double* pZetaA = new double[numOfContractionA];
            for (int i = 0; i < numOfContractionA; ++i) {
                pCoefA[i] = this->pOrbitalInfo_->getCoefficient(orbA, i);
                pZetaA[i] = this->pOrbitalInfo_->getExponent(orbA, i);
            }
            const TlPosition posA = this->pOrbitalInfo_->getPosition(orbA);
            double A[3];
            A[0] = posA.x();
            A[1] = posA.y();
            A[2] = posA.z();

            const int nqB = this->pOrbitalInfo_->getShellType(orbB);
            const int iAngB  = nqB * 2 + 1;
            const int numOfContractionB = this->pOrbitalInfo_->getCgtoContraction(orbB);
            double* pCoefB = new double[numOfContractionB];
            double* pZetaB = new double[numOfContractionB];
            for (int i = 0; i < numOfContractionB; ++i) {
                pCoefB[i] = this->pOrbitalInfo_->getCoefficient(orbB, i);
                pZetaB[i] = this->pOrbitalInfo_->getExponent(orbB, i);
            }
            const TlPosition posB = this->pOrbitalInfo_->getPosition(orbB);
            double B[3];
            B[0] = posB.x();
            B[1] = posB.y();
            B[2] = posB.z();

            for (std::size_t orbC = 0; orbC < numOfAuxDens;) {
                const int nqC = this->pOrbitalInfo_Density_->getShellType(orbC);
                const int iAngC  = nqC * 2 + 1;
                const double coefC = this->pOrbitalInfo_Density_->getCoefficient(orbC, 0);
                const double zetaC = this->pOrbitalInfo_Density_->getExponent(orbC, 0);
                const TlPosition posC = this->pOrbitalInfo_Density_->getPosition(orbC);
                double C[3];
                C[0] = posC.x();
                C[1] = posC.y();
                C[2] = posC.z();

                this->ericalc(nqA, nqB, nqC,
                              numOfContractionA, pCoefA, pZetaA, A,
                              numOfContractionB, pCoefB, pZetaB, B,
                              coefC, zetaC, C, ERI);

#pragma omp critical (DfEri_ericalcDT)
                {
                    int iww = 0;
                    for (int iAng = 0; iAng < iAngA; ++iAng) {
                        const std::size_t I = orbA + iAng;

                        for (int jAng = 0; jAng < iAngB; ++jAng) {
                            const std::size_t J = orbB + jAng;

                            if ((orbA != orbB) || (I >= J)) {
                                const double coef_P_pq = 2.0 * deltaPpq->get(I, J);
                                for (int kAng = 0; kAng < iAngC; ++kAng) {
                                    if (std::fabs(ERI[iww]) >= dCutValue) {
                                        const std::size_t K = orbC + kAng;

                                        (*pDeltaT)[K] += ERI[iww] * coef_P_pq;
                                    }
                                    ++iww;
                                }
                            } else {
                                iww += iAngC;
                            }
                        }
                    }
                }
                orbC += iAngC;
            }

            delete[] pCoefA;
            pCoefA = NULL;
            delete[] pCoefB;
            pCoefB = NULL;
            delete[] pZetaA;
            pZetaA = NULL;
            delete[] pZetaB;
            pZetaB = NULL;
        }
    }
}


void DfEri::ericalcDH(const std::vector<IJShellPair>& IJShellList, const TlVector& deltaRho, TlMatrixObject* deltaH)
{
    const std::size_t numOfAuxDens = this->m_nNumOfAux;

    const int maxNumOfPairs = IJShellList.size();
#pragma omp parallel
    {
        double ERI[(2 *MAX_SHELL+1) *(2 *MAX_SHELL+1) *(2 *MAX_SHELL+1)];
#pragma omp for
        for (int pairIndex = 0; pairIndex < maxNumOfPairs; ++pairIndex) {
            const std::size_t orbA = IJShellList[pairIndex].nIShell;
            const std::size_t orbB = IJShellList[pairIndex].nJShell;

            const int nqA = this->pOrbitalInfo_->getShellType(orbA);
            const int iAngA = nqA * 2 + 1;
            const int numOfContractionA = this->pOrbitalInfo_->getCgtoContraction(orbA);
            double* pCoefA = new double[numOfContractionA];
            double* pZetaA = new double[numOfContractionA];
            for (int i = 0; i < numOfContractionA; ++i) {
                pCoefA[i] = this->pOrbitalInfo_->getCoefficient(orbA, i);
                pZetaA[i] = this->pOrbitalInfo_->getExponent(orbA, i);
            }
            const TlPosition posA = this->pOrbitalInfo_->getPosition(orbA);
            double A[3];
            A[0] = posA.x();
            A[1] = posA.y();
            A[2] = posA.z();

            const int nqB = this->pOrbitalInfo_->getShellType(orbB);
            const int iAngB = nqB * 2 + 1;
            const int numOfContractionB = this->pOrbitalInfo_->getCgtoContraction(orbB);
            double* pCoefB = new double[numOfContractionB];
            double* pZetaB = new double[numOfContractionB];
            for (int i = 0; i < numOfContractionB; ++i) {
                pCoefB[i] = this->pOrbitalInfo_->getCoefficient(orbB, i);
                pZetaB[i] = this->pOrbitalInfo_->getExponent(orbB, i);
            }
            const TlPosition posB = this->pOrbitalInfo_->getPosition(orbB);
            double B[3];
            B[0] = posB.x();
            B[1] = posB.y();
            B[2] = posB.z();

            for (std::size_t orbC = 0; orbC < numOfAuxDens;) {
                const int nqC = this->pOrbitalInfo_Density_->getShellType(orbC);
                const int iAngC = nqC * 2 + 1;
                const double coefC = this->pOrbitalInfo_Density_->getCoefficient(orbC, 0);
                const double zetaC = this->pOrbitalInfo_Density_->getExponent(orbC, 0);
                const TlPosition posC = this->pOrbitalInfo_Density_->getPosition(orbC);
                double C[3];
                C[0] = posC.x();
                C[1] = posC.y();
                C[2] = posC.z();

                this->ericalc(nqA, nqB, nqC,
                              numOfContractionA, pCoefA, pZetaA, A,
                              numOfContractionB, pCoefB, pZetaB, B,
                              coefC, zetaC, C, ERI);

#pragma omp critical (DfEri_ericalcDH)
                {
                    int iww = 0;
                    for (int iAng = 0; iAng < iAngA; ++iAng) {
                        const std::size_t I = orbA + iAng;

                        for (int jAng = 0; jAng < iAngB; ++jAng) {
                            const std::size_t J = orbB + jAng;

                            if ((orbA != orbB) || (I >= J)) {
                                double value = 0.0;
                                for (int kAng = 0; kAng < iAngC; ++kAng) {
                                    const std::size_t K = orbC + kAng;

                                    value += ERI[iww] * deltaRho[K];
                                    ++iww;
                                }
                                deltaH->add(I, J, value);
                            } else {
                                iww += iAngC;
                            }
                        }
                    }
                }
                orbC += iAngC;
            }

            delete[] pCoefA;
            pCoefA = NULL;
            delete[] pCoefB;
            pCoefB = NULL;
            delete[] pZetaA;
            pZetaA = NULL;
            delete[] pZetaB;
            pZetaB = NULL;
        }
    }
}


bool DfEri::isCutoffUsingSchwartzInequality(const TlOrbitalInfo& orbitalInfo,
                                            const std::size_t p, const std::size_t q,
                                            const double dCutoffValue) const
{
    bool bAnswer = true;

    // SS OVERLAP
    const TlPosition A = orbitalInfo.getPosition(p);
    const TlPosition B = orbitalInfo.getPosition(q);
    const double absQ = A.squareDistanceFrom(B);
    const double CPAI = M_PI * sqrt(M_PI);
    double STS = 0.0;

    const int numOfContractionA = orbitalInfo.getCgtoContraction(p);
    const int numOfContractionB = orbitalInfo.getCgtoContraction(q);
    for (int ipA = 0; ipA < numOfContractionA; ++ipA) {
        const double ZAW = orbitalInfo.getExponent(p, ipA);
        const double coeff = CPAI * orbitalInfo.getCoefficient(p, ipA);

        for (int ipB = 0; ipB < numOfContractionB; ++ipB) {
            const double ZBW = orbitalInfo.getExponent(q, ipB);
            const double ZP = ZAW + ZBW;
            const double ZPIW = 1.0/ZP;
            STS += coeff * ZPIW * sqrt(ZPIW) * exp(-absQ * ZAW * ZBW * ZPIW) * orbitalInfo.getCoefficient(q, ipB);
        }
    }
    // SS OVERLAP END

    if (std::fabs(STS) >= dCutoffValue) {
        bAnswer = false;
    }

    return bAnswer;
}

/*****************************************************************************/
/*   ERI calculation for [pq|Alpha] ( [AB|C] )                   */
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
/*  FPAI = 2 * PAI**(5/2)                            */
/*****************************************************************************/
void DfEri::ericalc(const int nqA, const int nqB, const int nqC,
                    const int npA, const double* CA, const double* ZA, const double* A,
                    const int npB, const double* CB, const double* ZB, const double* B,
                    const double CC, const double ZC, const double* C,
                    double* ERI)
{
    assert((0 <= nqA) && (nqA <= 3));
    assert((0 <= nqB) && (nqB <= 3));
    assert((0 <= nqC) && (nqC <= 3));
    assert(nqA >= nqB);
    assert((npA <= 10) && (npB <= 10));

    const int np    = npA * npB;
    const int iAngA = nqA * 2 + 1;
    const int iAngB = nqB * 2 + 1;
    const int iAngC = nqC * 2 + 1;
    const int nAng  = iAngA * iAngB * iAngC;

    double absQ = 0.0;
    for (int i = 0; i < 3; ++i) {
        const double ABW = A[i] - B[i];
        absQ += ABW * ABW;
    }

    const int contractions = npA * npB;
    std::vector<TlPosition> P(contractions);
    std::vector<TlPosition> PA(contractions);
    std::vector<TlPosition> PB(contractions);

    double* Px = new double[contractions];
    double* Py = new double[contractions];
    double* Pz = new double[contractions];
    double* PAx = new double[contractions];
    double* PAy = new double[contractions];
    double* PAz = new double[contractions];
    double* PBx = new double[contractions];
    double* PBy = new double[contractions];
    double* PBz = new double[contractions];

    const double Gammainv = 1.0 / ZC;
    double* HP = new double[contractions];
    double* Zpi = new double[contractions];

    double** ERP = new double*[contractions];
    for (int i = 0; i < contractions; ++i) {
        ERP[i] = new double[(MAX_SHELL * 2 +1) *(MAX_SHELL * 2 +1) *(MAX_SHELL * 2 +1)];
    }

    {
        std::size_t ip = 0;
        for (int ipA = 0; ipA < npA; ++ipA) {
            const double ZAW  = ZA[ipA];

            for (int ipB = 0; ipB < npB; ++ipB) {
                const double ZBW  = ZB[ipB];
                const double ZP   = ZAW + ZBW;
                const double ZPIW = 1.0 / ZP;
                Zpi[ip] = ZPIW;
                HP[ip]  = ZPIW * sqrt(ZPIW) * exp(-absQ * ZAW * ZBW * ZPIW) * CA[ipA] * CB[ipB];

                const double px = (ZAW * A[0] + ZBW * B[0])*ZPIW;
                const double py = (ZAW * A[1] + ZBW * B[1])*ZPIW;
                const double pz = (ZAW * A[2] + ZBW * B[2])*ZPIW;
                Px[ip] = px;
                Py[ip] = py;
                Pz[ip] = pz;
                PAx[ip] = px - A[0];
                PAy[ip] = py - A[1];
                PAz[ip] = pz - A[2];
                PBx[ip] = px - B[0];
                PBy[ip] = py - B[1];
                PBz[ip] = pz - B[2];

                ++ip;
            }
        }
    }

    const double Cx = C[0];
    const double Cy = C[1];
    const double Cz = C[2];

    const int type = nqA * 16 + nqB * 4 + nqC;
    switch (type) {
    case 0:
        this->eriSSS(np, Px, Py, Pz, Cx, Cy, Cz, CC, Gammainv, Zpi, HP, ERP);
        break;

    case 1:
        this->eriSSP(np, Px, Py, Pz, Cx, Cy, Cz, CC, Gammainv, Zpi, HP, ERP);
        break;

    case 2:
        this->eriSSD(np, Px, Py, Pz, Cx, Cy, Cz, CC, Gammainv, Zpi, HP, ERP);
        break;

    case 16:
        this->eriPSS(np, Px, Py, Pz, PAx, PAy, PAz, Cx, Cy, Cz, CC, Gammainv, Zpi, HP, ERP);
        break;

    case 17:
        this->eriPSP(np, Px, Py, Pz, PAx, PAy, PAz, Cx, Cy, Cz, CC, Gammainv, Zpi, HP, ERP);
        break;

    case 18:
        this->eriPSD(np, Px, Py, Pz, PAx, PAy, PAz, Cx, Cy, Cz, CC, Gammainv, Zpi, HP, ERP);
        break;

    case 20:
        this->eriPPS(np, Px, Py, Pz, PAx, PAy, PAz, PBx, PBy, PBz, Cx, Cy, Cz, CC, Gammainv, Zpi, HP, ERP);
        break;

    case 21:
        this->eriPPP(np, Px, Py, Pz, PAx, PAy, PAz, PBx, PBy, PBz, Cx, Cy, Cz, CC, Gammainv, Zpi, HP, ERP);
        break;

    case 22:
        this->eriPPD(np, Px, Py, Pz, PAx, PAy, PAz, PBx, PBy, PBz, Cx, Cy, Cz, CC, Gammainv, Zpi, HP, ERP);
        break;

    case 32:
        this->eriDSS(np, Px, Py, Pz, PAx, PAy, PAz, Cx, Cy, Cz, CC, Gammainv, Zpi, HP, ERP);
        break;

    case 33:
        this->eriDSP(np, Px, Py, Pz, PAx, PAy, PAz, Cx, Cy, Cz, CC, Gammainv, Zpi, HP, ERP);
        break;

    case 34:
        this->eriDSD(np, Px, Py, Pz, PAx, PAy, PAz, Cx, Cy, Cz, CC, Gammainv, Zpi, HP, ERP);
        break;

    case 36:
        this->eriDPS(np, Px, Py, Pz, PAx, PAy, PAz, PBx, PBy, PBz, Cx, Cy, Cz, CC, Gammainv, Zpi, HP, ERP);
        break;

    case 37:
        this->eriDPP(np, Px, Py, Pz, PAx, PAy, PAz, PBx, PBy, PBz, Cx, Cy, Cz, CC, Gammainv, Zpi, HP, ERP);
        break;

    case 38:
        this->eriDPD(np, Px, Py, Pz, PAx, PAy, PAz, PBx, PBy, PBz, Cx, Cy, Cz, CC, Gammainv, Zpi, HP, ERP);
        break;

    case 40:
        this->eriDDS(np, Px, Py, Pz, PAx, PAy, PAz, PBx, PBy, PBz, Cx, Cy, Cz, CC, Gammainv, Zpi, HP, ERP);
        break;

    case 41:
        this->eriDDP(np, Px, Py, Pz, PAx, PAy, PAz, PBx, PBy, PBz, Cx, Cy, Cz, CC, Gammainv, Zpi, HP, ERP);
        break;

    case 42:
        this->eriDDD(np, Px, Py, Pz, PAx, PAy, PAz, PBx, PBy, PBz, Cx, Cy, Cz, CC, Gammainv, Zpi, HP, ERP);
        break;

    default:
        abort();
        break;
    }

    // Contraction loop
//   for (int i = 0; i < nAng; ++i){
//     ERI[i] = 0.0;
//   }

    for (int i = 0; i < nAng; ++i) {
        double tmp = 0.0;

        for (int ip = 0; ip < np; ++ip) {
            tmp += ERP[ip][i];
        }

        ERI[i] = tmp;
    }

    for (int i = 0; i < contractions; ++i) {
        delete[] ERP[i];
    }
    delete[] ERP;

    delete[] HP;
    delete[] Zpi;
    delete[] Px;
    delete[] Py;
    delete[] Pz;
    delete[] PAx;
    delete[] PAy;
    delete[] PAz;
    delete[] PBx;
    delete[] PBy;
    delete[] PBz;
}


/*****************************************************************************/
/*   make table for ERI [pq|Alpha]                       */
/*  natom  = total number of atoms                       */
/*  ncgto  =         orbital CGTOs                   */
/*  ndgto  =         density function GTOs               */
/*  norbf  =         orbital basis functions             */
/*  ndensf =         density functions               */
/*****************************************************************************/
void DfEri::makeTable()
{
    const std::size_t numOfAOs = this->m_nNumOfAOs;
    this->shellList_.clear();
    this->shellList_.resize(3); // s, p, d support
    {
        std::size_t shell = 0;
        while (shell < numOfAOs) {
            const int shellType = this->pOrbitalInfo_->getShellType(shell);
            this->shellList_[shellType].push_back(shell);
            shell += shellType * 2 + 1;
        }
    }

    const index_type numOfAuxDens = this->m_nNumOfAux;
    assert(this->pOrbitalInfo_Density_->getNumOfOrbitals() == numOfAuxDens);
    this->shellList_Dens_.clear();
    this->shellList_Dens_.resize(3); // s, p, d support
    {
        index_type shell = 0;
        while (shell < numOfAuxDens) {
            const int shellType = this->pOrbitalInfo_Density_->getShellType(shell);
            this->shellList_Dens_[shellType].push_back(shell);
            shell += shellType * 2 + 1;
        }
    }
}


/*****************************************************************************/
/*   aux set for ERI [pq|Alpha]                          */
/*  FPAI   = 2*PAI**2.5                          */
/*  D33    = 1/0.03                              */
/*  SQR3I  = 1/sqrt(3)                           */
/*  SQR3I2 = 2*SQR3I                             */
/*****************************************************************************/
int DfEri::auxSet()
{
    //int l;
    double TF0[25], D[6];
    double DD, GA5;
    int nk[6];

    // TF Data
    TF0[ 0] = 33.0;
    TF0[ 1] = 37.0;
    TF0[ 2] = 41.0;
    TF0[ 3] = 43.0;
    TF0[ 4] = 46.0;
    TF0[ 5] = 49.0;
    TF0[ 6] = 51.0;
    TF0[ 7] = 54.0;
    TF0[ 8] = 56.0;
    TF0[ 9] = 58.0;
    TF0[10] = 61.0;
    TF0[11] = 63.0;
    TF0[12] = 66.0;
    TF0[13] = 68.0;
    TF0[14] = 70.0;
    TF0[15] = 72.0;
    TF0[16] = 74.0;
    TF0[17] = 80.0;
    TF0[18] = 85.0;
    TF0[19] = 90.0;
    TF0[20] = 95.0;
    TF0[21] = 95.0;
    TF0[22] = 98.0;
    TF0[23] = 99.0;
    TF0[24] = 99.0;

    nk[0] =  -1;
    nk[2] = 18;
    nk[4] = -48;
    nk[1] = nk[3] = nk[5] = 0;

    DD    = 0.015;
    D[1]  = DD*DD;
    D[3]  = D[1]*D[1];
    D[5]  = D[1]*D[3];

    for (int k = 0; k < 25; ++k) {
        TF[k] = TF0[k];
    }

    for (int k = 1; k < 25; ++k) {
        RMI[k-1] = 1.0 / (double)(2*k+1);
    }

    GA[0] = 1.0;
    for (int k = 1; k < 6; ++k) {
        GA[k] = GA[k-1]/(double) k;
    }

    GA5  = GA[5] / 192.0;             // 192 = 6*32
    for (int k = 0; k < 6; ++k) {
        if (nk[k] != 0) {
            GA[k] -= nk[k]*GA5*D[5-k];
        }
    }

    // The following loop may be meaningless, because the same loop exists
    // in member function, fmtSet.
    for (int k = 0; k < 1901; ++k) {
        EDAT[k] = exp(-0.03*k);
    }

    this->fmtSet();

    // copy *FDAT[6]  to *SDAT
    for (int k = 0; k < 1901; ++k) {
        SDAT[k] = FDAT[6][k];
    }

    this->fmtRecursive(14, 1901, 5);

    //
    int mp  = 2;
    int mnp = 1;
    for (int m = 8; m >= 0; m--) {
        int mm = m*6;
        double gm = 1.0;

        if (m != 8) {
            mp--;
            if (mp == 0) {
                mp  = 7;
            }

            mnp--;
            if (mnp == 0) {
                mnp = 7;
            }

            // fmtRecursive(m+1, 1901, 0);
            const double DMI = 1.0 / (2.0 * m + 1.0);
            for (int k = 0; k < 1901; ++k) {
                const double t2 = 0.06*(double) k;
                FDAT[mnp-1][k] = DMI*(t2*FDAT[mp-1][k] + EDAT[k]);
            }
            // fmtRecursive end
        }

        int jm = mnp - 1;
        if (jm == 0) {
            jm = 7;
        }
        for (int  k = 0; k < 6; ++k) {
            int mk = mm + k;
            if (k != 0) {
                gm /= (double) k;
            }
            int nkk   = nk[k];
            double dk = D[5-k];
            int jk    = mnp + k;
            int jj    = jk / 8;   // check here
            jk       -= jj*7;     // check here
            if (nkk == 0) {
                for (int l = 0; l < 1901; ++l) {
                    ADAT[mk][l] = gm*FDAT[jk-1][l];
                }
                continue;
            }
            const double gdk = nkk*dk*GA5;

            for (int l = 0; l < 1901; ++l) {
                ADAT[mk][l] = gm*FDAT[jk-1][l] - gdk*FDAT[jm-1][l];
            }
        }
    }

    return 0;

}

/*****************************************************************************/
/*  FMT calculation table set (m = 14)                   */
/*  Should calculate "int double" (real*16).                 */
/*****************************************************************************/
void DfEri::fmtSet()
{
    int n1 =0;
    int n2 =0;
    int n3 =0;
    double T[1901], T2[1901];
    for (int k = 0; k < 1901; ++k) {
        T[k]  = 0.03*(double)k;
        T2[k] = 2*T[k];
        if (T[k] <=  2.0) {
            n1 = k;
        }
        if (T[k] <= 10.0) {
            n2 = k;
        }
        if (T[k] <= 20.0) {
            n3 = k;
        }
    }

    double FMT[1901], Phai[1901];

    for (int k = 0; k < 1901; ++k) {
        this->EDAT[k] = exp(-T[k]);
        FMT[k]  = 0.0;
        Phai[k] = 0.0;
    }

    for (int m = 50; m >= 14; m--) {
        const double DMI = 1.0 / (2.0*m + 1.0);
        for (int k = 0; k < n1; ++k)
            FMT[k] = DMI*(T2[k]*FMT[k] + EDAT[k]);
    }

    for (int m = 70; m >= 14; m--) {
        const double DMI = 1.0 / (2.0*m + 1.0);
        for (int k = n1; k < n2; ++k)
            FMT[k] = DMI*(T2[k]*FMT[k] + EDAT[k]);
    }

    for (int m = 120; m >= 14; m--) {
        const double DMI = 1.0 / (2.0*m + 1.0);
        for (int k = n2; k < n3; ++k) {
            FMT[k] = DMI*(T2[k]*FMT[k] + EDAT[k]);
        }
    }

    for (int k = n3; k < 1901; ++k) {
        FMT[k] = 0.5*sqrt(M_PI /T[k]);
    }

    for (int m = 1; m <= 14; ++m) {
        for (int k = n3; k < 1901; ++k) {
            FMT[k] *= (2.0*m - 1.0) / T2[k];
        }
    }

    for (int k = n3; k < 1901; ++k) {
        int ms = (int)(T[k] + 0.5);
        ms *= -1;
        Phai[k] = EDAT[k] / (T2[k] - 2.0*ms);
        for (int m = ms+1; m <= 14; ++m) {
            Phai[k] = ((2.0*m - 1.0)*Phai[k] + EDAT[k]) / T2[k];
        }
    }

    for (int k = n3; k < 1901; ++k) {
        FMT[k] -= Phai[k];
    }

    for (int k = 0; k < 1901; ++k) {
        this->FDAT[6][k] = FMT[k];      // 0-5 are work area.
    }

    return;

}

/*****************************************************************************/
/*  FMT backward recursive relation                      */
/*****************************************************************************/
void DfEri::fmtRecursive(int m, int n, int mm)
{
    assert(mm >= 1);
    {
        const double DMI = 1.0 / (2.0*m - 1.0);
        for (int k = 0; k < n; ++k) {
            const double t2 = 0.06*(double) k;
            this->FDAT[mm][k] = DMI * (t2*SDAT[k] + EDAT[k]);
        }
    }

    int kk = m - 1;
    for (int k = 1; k <= mm; ++k) {
        kk--;
        const double DMI = 1.0 / (2.0*kk + 1.0);
        for (int l = 0; l < n; ++l) {
            const double t2 = 0.06*(double) l;
            this->FDAT[mm-k][l] = DMI*(t2*FDAT[mm-k+1][l] + EDAT[l]);
        }
    }

    return;
}

/*****************************************************************************/
/*  Primitive 3-center ERI evaluation [ AB | C ] for type [ SS | S ]         */
/*  np      ; number of first primitive pair data            */
/*  Zpi     ; 1/(Alpha+Beta) from orbital exponents Alpha, Beta  */
/*  px, py, pz  ; (Alpha*ax * Beta*bx) / (Alpha+Beta) , ...          */
/*  HP      ; Zeta**(-3/2) * exp(-(Alpha*Beta/Zeta) * (A-B)**2   */
/*  Gammainv    ; inverse of orbital exponent Gamma for last GTO     */
/*  cx, cy, cz  ; coordinates of last GTO C              */
/*  cc      ; coefficient of last GTO C              */
/*  ERP     ; output primitive ERI                   */
/*  FPAI = 2 * PAI**(5/2)                            */
/*****************************************************************************/
void DfEri::eriSSS(const int np,
                   const double* px, const double* py, const double* pz,
                   const double cx, const double cy, const double cz,
                   const double cc, const double Gammainv, const double* Zpi, const double* HP, double** ERP)
{
    const double tfw = TF[0];
    const double hc  = FPAI*Gammainv*sqrt(Gammainv)*cc;

    for (int i = 0; i < np; ++i) {
        const double Roo = 1.0 / (Zpi[i] + Gammainv);
        const double pcx = px[i] - cx;
        const double pcy = py[i] - cy;
        const double pcz = pz[i] - cz;
        const double t = Roo * (pcx*pcx + pcy*pcy + pcz*pcz);

        double f0;
        if (t <= tfw) {
            const int it = (int)((t+0.015) * D33);
            const double dt = 0.03*it - t;
            f0 = ((((ADAT[5][it]  *dt + ADAT[4][it]) * dt
                    + ADAT[3][it])*dt + ADAT[2][it]) * dt
                  + ADAT[1][it])*dt + ADAT[0][it];
        } else {
            const double tinv = 1.0/t;
            f0 = sqrt(M_PI * tinv) * 0.5;
        }

        ERP[i][0] = hc * sqrt(Roo) * HP[i] * f0;
    }
}

/*****************************************************************************/
/*  Primitive 3-center ERI evaluation [ AB | C ] for type [ PS | S ]         */
/*  np      ; number of first primitive pair data            */
/*  Zpi     ; 1/(Alpha+Beta) from orbital exponents Alpha, Beta  */
/*  px, py, pz  ; (Alpha*ax * Beta*bx) / (Alpha+Beta) , ...          */
/*  pax, pay, paz   ; px-ax, py-ay, pz-az                    */
/*  HP      ; Zeta**(-3/2) * exp(-(Alpha*Beta/Zeta) * (A-B)**2   */
/*  Gammainv    ; inverse of orbital exponent Gamma for last GTO     */
/*  cx, cy, cz  ; coordinates of last GTO C              */
/*  cc      ; coefficient of last GTO C              */
/*  ERP     ; output primitive ERI                   */
/*              ERP[0,1,2,...,np-1][0] for [ PxS | S ]       */
/*              ERP[0,1,2,...,np-1][1] for [ PyS | S ]       */
/*              ERP[0,1,2,...,np-1][2] for [ PzS | S ]       */
/*  FPAI = 2 * PAI**(5/2)                            */
/*****************************************************************************/
void DfEri::eriPSS(const int np,
                   const double* px, const double* py, const double* pz,
                   const double* pax, const double* pay, const double* paz,
                   const double cx, const double cy, const double cz,
                   const double cc, const double Gammainv, const double* Zpi, const double* HP, double** ERP)
{
    const double tfw = TF[1];
    const double hc = FPAI*Gammainv*sqrt(Gammainv)*cc;

    for (int i = 0; i < np; ++i) {
        const double Roo = 1.0 / (Zpi[i] + Gammainv);
        const double pcx = px[i] - cx;
        const double pcy = py[i] - cy;
        const double pcz = pz[i] - cz;
        const double t   = Roo * (pcx*pcx + pcy*pcy + pcz*pcz);

        double f0, f1;
        if (t <= tfw) {
            const int it  = (int)((t+0.015) * D33);
            const double dt  = 0.03*it - t;
            f1  = ((((ADAT[11][it]  *dt + ADAT[10][it]) * dt
                     + ADAT[ 9][it])*dt + ADAT[ 8][it]) * dt
                   + ADAT[ 7][it])*dt + ADAT[ 6][it];
            const double eed = ((((GA[5]         *dt + GA[4]) * dt
                                  + GA[3])*dt + GA[2]) * dt
                                +   GA[1])*dt + GA[0];
            const double ee  = EDAT[it] * eed;
            f0  = 2.0*t*f1 + ee;
        } else {
            const double tinv = 1.0/t;
            f0   = sqrt(M_PI * tinv)*0.5;
            f1   = 0.5*f0*tinv;
        }

        const double hpq   = hc*sqrt(Roo)*HP[i];
        const double sss0  = hpq*f0;
        const double Roof1 = Zpi[i]*Roo*hpq*f1;

        ERP[i][0] = pax[i]*sss0 - pcx*Roof1;
        ERP[i][1] = pay[i]*sss0 - pcy*Roof1;
        ERP[i][2] = paz[i]*sss0 - pcz*Roof1;
    }
}


/*****************************************************************************/
/*                                       */
/*  Primitive 3-center ERI evaluation [ AB | C ] for type [ SS | P ]         */
/*                                       */
/*  np      ; number of first primitive pair data            */
/*  Zpi     ; 1/(Alpha+Beta) from orbital exponents Alpha, Beta  */
/*  px, py, pz  ; (Alpha*ax * Beta*bx) / (Alpha+Beta) , ...          */
/*  HP      ; Zeta**(-3/2) * exp(-(Alpha*Beta/Zeta) * (A-B)**2   */
/*  Gammainv    ; inverse of orbital exponent Gamma for last GTO     */
/*  cx, cy, cz  ; coordinates of last GTO C              */
/*  cc      ; coefficient of last GTO C              */
/*  ERP     ; output primitive ERI                   */
/*              ERP[0,1,2,...,np-1][0] for [ SS | Px ]       */
/*              ERP[0,1,2,...,np-1][1] for [ SS | Py ]       */
/*              ERP[0,1,2,...,np-1][2] for [ SS | Pz ]       */
/*                                       */
/*  FPAI = 2 * PAI**(5/2)                            */
/*                                       */
/*****************************************************************************/
void DfEri::eriSSP(const int np,
                   const double* px, const double* py, const double* pz,
                   const double cx, const double cy, const double cz,
                   const double cc, const double Gammainv, const double* Zpi, const double* HP, double** ERP)
{
    const double tfw = TF[1];
    const double hc  = FPAI*Gammainv*sqrt(Gammainv)*cc;

    for (int i = 0; i < np; ++i) {
        const double Roo = 1.0 / (Zpi[i] + Gammainv);
        const double pcx = px[i] - cx;
        const double pcy = py[i] - cy;
        const double pcz = pz[i] - cz;
        const double t = Roo * (pcx*pcx + pcy*pcy + pcz*pcz);

        double f1;
        if (t <= tfw) {
            const int it  = (int)((t+0.015) * D33);
            const double dt  = 0.03*it - t;
            f1  = ((((ADAT[11][it]  *dt + ADAT[10][it]) * dt
                     + ADAT[ 9][it])*dt + ADAT[ 8][it]) * dt
                   + ADAT[ 7][it])*dt + ADAT[ 6][it];
        } else {
            const double tinv = 1.0/t;
            f1   = 0.25*tinv*sqrt(M_PI * tinv);
        }

        const double Roof1 = Roo*Gammainv*f1*hc*sqrt(Roo)*HP[i];

        ERP[i][0] = pcx*Roof1;
        ERP[i][1] = pcy*Roof1;
        ERP[i][2] = pcz*Roof1;
    }
}


/*****************************************************************************/
/*                                       */
/*  Primitive 3-center ERI evaluation [ AB | C ] for type [ PP | S ]         */
/*                                       */
/*  np      ; number of first primitive pair data            */
/*  Zpi     ; 1/(Alpha+Beta) from orbital exponents Alpha, Beta  */
/*  px, py, pz  ; (Alpha*ax * Beta*bx) / (Alpha+Beta) , ...          */
/*  pax, pay, paz   ; px-ax, py-ay, pz-az                    */
/*  pbx, pby, pbz   ; px-bx, py-by, pz-bz                    */
/*  HP      ; Zeta**(-3/2) * exp(-(Alpha*Beta/Zeta) * (A-B)**2   */
/*  Gammainv    ; inverse of orbital exponent Gamma for last GTO     */
/*  cx, cy, cz  ; coordinates of last GTO C              */
/*  cc      ; coefficient of last GTO C              */
/*  ERP     ; output primitive ERI                   */
/*              ERP[0,1,2,...,np-1][0] for [ PxPx | S ]      */
/*              ERP[0,1,2,...,np-1][1] for [ PxPy | S ]      */
/*              ERP[0,1,2,...,np-1][2] for [ PxPz | S ]      */
/*              ERP[0,1,2,...,np-1][3] for [ PyPx | S ]      */
/*              ERP[0,1,2,...,np-1][4] for [ PyPy | S ]      */
/*              ERP[0,1,2,...,np-1][5] for [ PyPz | S ]      */
/*              ERP[0,1,2,...,np-1][6] for [ PzPx | S ]      */
/*              ERP[0,1,2,...,np-1][7] for [ PzPy | S ]      */
/*              ERP[0,1,2,...,np-1][8] for [ PzPz | S ]      */
/*                                       */
/*  FPAI = 2 * PAI**(5/2)                            */
/*                                       */
/*****************************************************************************/
void DfEri::eriPPS(const int np,
                   const double* px, const double* py, const double* pz,
                   const double* pax, const double* pay, const double* paz,
                   const double* pbx, const double* pby, const double* pbz,
                   const double cx, const double cy, const double cz,
                   const double cc, const double Gammainv, const double* Zpi, const double* HP, double** ERP)
{
    const double tfw  = TF[2];
    const double rmi0 = RMI[0];
    const double hc   = FPAI*Gammainv*sqrt(Gammainv)*cc;

    for (int i = 0; i < np; ++i) {
        const double Roo = 1.0 / (Zpi[i] + Gammainv);
        const double pcx = px[i] - cx;
        const double pcy = py[i] - cy;
        const double pcz = pz[i] - cz;
        const double t = Roo * (pcx*pcx + pcy*pcy + pcz*pcz);

        double f0, f1, f2;
        if (t <= tfw) {
            const int it  = (int)((t+0.015) * D33);
            const double dt  = 0.03*it - t;
            f2  = ((((ADAT[17][it]  *dt + ADAT[16][it]) * dt
                     + ADAT[15][it])*dt + ADAT[14][it]) * dt
                   + ADAT[13][it])*dt + ADAT[12][it];
            const double eed = ((((GA[5]         *dt + GA[4]) * dt
                                  + GA[3])*dt + GA[2]) * dt
                                +   GA[1])*dt + GA[0];
            const double ee  = EDAT[it] * eed;
            const double t2  = 2.0*t;
            f1  = (t2*f2 + ee)*rmi0;
            f0  = t2*f1 + ee;
        } else {
            const double tinv = 1.0/t;
            f0   = sqrt(M_PI * tinv)*0.5;
            f1   = 0.5*f0*tinv;
            f2   = 1.5*f1*tinv;
        }

        const double hpq   = hc*sqrt(Roo)*HP[i];
        const double sss0  = hpq*f0;
        const double sss1  = hpq*f1;
        const double sss2  = hpq*f2;
        // [ PS | S ]
        const double Rooz  = Zpi[i]*Roo;
        const double pgx   = Rooz*pcx;
        const double pgy   = Rooz*pcy;
        const double pgz   = Rooz*pcz;
        const double xss0  = pax[i]*sss0 - pgx*sss1;
        const double yss0  = pay[i]*sss0 - pgy*sss1;
        const double zss0  = paz[i]*sss0 - pgz*sss1;
        const double xss1  = pax[i]*sss1 - pgx*sss2;
        const double yss1  = pay[i]*sss1 - pgy*sss2;
        const double zss1  = paz[i]*sss1 - pgz*sss2;
        // [ PP | S ]
        const double Zp5   = 0.5*Zpi[i];
        const double Zps01 = Zp5*(sss0 - Rooz*sss1);

        ERP[i][0] = pbx[i]*xss0 - pgx*xss1 + Zps01;
        ERP[i][1] = pby[i]*xss0 - pgy*xss1;
        ERP[i][2] = pbz[i]*xss0 - pgz*xss1;
        ERP[i][3] = pbx[i]*yss0 - pgx*yss1;
        ERP[i][4] = pby[i]*yss0 - pgy*yss1 + Zps01;
        ERP[i][5] = pbz[i]*yss0 - pgz*yss1;
        ERP[i][6] = pbx[i]*zss0 - pgx*zss1;
        ERP[i][7] = pby[i]*zss0 - pgy*zss1;
        ERP[i][8] = pbz[i]*zss0 - pgz*zss1 + Zps01;
    }
}


/*****************************************************************************/
/*                                       */
/*  Primitive 3-center ERI evaluation [ AB | C ] for type [ PS | P ]         */
/*                                       */
/*  np      ; number of first primitive pair data            */
/*  Zpi     ; 1/(Alpha+Beta) from orbital exponents Alpha, Beta  */
/*  px, py, pz  ; (Alpha*ax * Beta*bx) / (Alpha+Beta) , ...          */
/*  pax, pay, paz   ; px-ax, py-ay, pz-az                    */
/*  HP      ; Zeta**(-3/2) * exp(-(Alpha*Beta/Zeta) * (A-B)**2   */
/*  Gammainv    ; inverse of orbital exponent Gamma for last GTO     */
/*  cx, cy, cz  ; coordinates of last GTO C              */
/*  cc      ; coefficient of last GTO C              */
/*  ERP     ; output primitive ERI                   */
/*              ERP[0,1,2,...,np-1][0] for [ PxS | Px ]      */
/*              ERP[0,1,2,...,np-1][1] for [ PyS | Px ]      */
/*              ERP[0,1,2,...,np-1][2] for [ PzS | Px ]      */
/*              ERP[0,1,2,...,np-1][3] for [ PxS | Py ]      */
/*              ERP[0,1,2,...,np-1][4] for [ PyS | Py ]      */
/*              ERP[0,1,2,...,np-1][5] for [ PzS | Py ]      */
/*              ERP[0,1,2,...,np-1][6] for [ PxS | Pz ]      */
/*              ERP[0,1,2,...,np-1][7] for [ PyS | Pz ]      */
/*              ERP[0,1,2,...,np-1][8] for [ PzS | Pz ]      */
/*                                       */
/*  FPAI = 2 * PAI**(5/2)                            */
/*                                       */
/*****************************************************************************/
void DfEri::eriPSP(int np,
                   const double* px, const double* py, const double* pz,
                   const double* pax, const double* pay, const double* paz,
                   const double cx, const double cy, const double cz,
                   const double cc, const double Gammainv, const double* Zpi, const double* HP, double** ERP)
{
    const double tfw  = TF[2];
    const double rmi0 = RMI[0];
    const double hc = FPAI*Gammainv*sqrt(Gammainv)*cc;

    for (int i = 0; i < np; ++i) {
        const double Roo = 1.0 / (Zpi[i] + Gammainv);
        const double pcx = px[i] - cx;
        const double pcy = py[i] - cy;
        const double pcz = pz[i] - cz;
        const double t   = Roo * (pcx*pcx + pcy*pcy + pcz*pcz);

        double f1, f2;
        if (t <= tfw) {
            const int it  = (int)((t+0.015) * D33);
            const double dt  = 0.03*it - t;
            f2  = ((((ADAT[17][it]  *dt + ADAT[16][it]) * dt
                     + ADAT[15][it])*dt + ADAT[14][it]) * dt
                   + ADAT[13][it])*dt + ADAT[12][it];
            const double eed = ((((GA[5]         *dt + GA[4]) * dt
                                  + GA[3])*dt + GA[2]) * dt
                                +   GA[1])*dt + GA[0];
            const double ee  = EDAT[it] * eed;
            const double t2  = 2.0*t;
            f1  = (t2*f2 + ee)*rmi0;
        } else {
            const double tinv = 1.0/t;
            f1   = 0.25*tinv*sqrt(M_PI * tinv);
            f2   = 1.5*f1*tinv;
        }

        const double hpq   = hc*sqrt(Roo)*HP[i];
        const double sss1  = hpq*f1;
        const double sss2  = hpq*f2;
        // [ PS | S ]
        const double Roof2 = Zpi[i]*Roo*sss2;
        const double xss1  = pax[i]*sss1 - pcx*Roof2;
        const double yss1  = pay[i]*sss1 - pcy*Roof2;
        const double zss1  = paz[i]*sss1 - pcz*Roof2;

        const double RooG  = Roo*Gammainv;
        const double Gcx   = RooG*pcx;
        const double Gcy   = RooG*pcy;
        const double Gcz   = RooG*pcz;
        const double RGf1  = 0.5*Zpi[i]*RooG*sss1;
        // [ PS | P ]
        ERP[i][0] = Gcx*xss1 + RGf1;
        ERP[i][1] = Gcy*xss1;
        ERP[i][2] = Gcz*xss1;
        ERP[i][3] = Gcx*yss1;
        ERP[i][4] = Gcy*yss1 + RGf1;
        ERP[i][5] = Gcz*yss1;
        ERP[i][6] = Gcx*zss1;
        ERP[i][7] = Gcy*zss1;
        ERP[i][8] = Gcz*zss1 + RGf1;
    }
}


/*****************************************************************************/
/*                                       */
/*  Primitive 3-center ERI evaluation [ AB | C ] for type [ PP | P ]         */
/*                                       */
/*  np      ; number of first primitive pair data            */
/*  Zpi     ; 1/(Alpha+Beta) from orbital exponents Alpha, Beta  */
/*  px, py, pz  ; (Alpha*ax * Beta*bx) / (Alpha+Beta) , ...          */
/*  pax, pay, paz   ; px-ax, py-ay, pz-az                    */
/*  pbx, pby, pbz   ; px-bx, py-by, pz-bz                    */
/*  HP      ; Zeta**(-3/2) * exp(-(Alpha*Beta/Zeta) * (A-B)**2   */
/*  Gammainv    ; inverse of orbital exponent Gamma for last GTO     */
/*  cx, cy, cz  ; coordinates of last GTO C              */
/*  cc      ; coefficient of last GTO C              */
/*  ERP     ; output primitive ERI                   */
/*              ERP[0,1,2,...,np-1][0] for [ PxPx | Px ]     */
/*              ERP[0,1,2,...,np-1][1] for [ PxPy | Px ]     */
/*              ERP[0,1,2,...,np-1][2] for [ PxPz | Px ]     */
/*              ERP[0,1,2,...,np-1][3] for [ PyPx | Px ]     */
/*              ERP[0,1,2,...,np-1][4] for [ PyPy | Px ]     */
/*              ERP[0,1,2,...,np-1][5] for [ PyPz | Px ]     */
/*              ERP[0,1,2,...,np-1][6] for [ PzPx | Px ]     */
/*              ERP[0,1,2,...,np-1][7] for [ PzPy | Px ]     */
/*              ERP[0,1,2,...,np-1][8] for [ PzPz | Px ]     */
/*                  .............                */
/*                                       */
/*  FPAI = 2 * PAI**(5/2)                            */
/*                                       */
/*****************************************************************************/
void DfEri::eriPPP(const int np,
                   const double* px, const double* py, const double* pz,
                   const double* pax, const double* pay, const double* paz,
                   const double* pbx, const double* pby, const double* pbz,
                   const double cx, const double cy, const double cz,
                   const double cc, const double Gammainv, const double* Zpi, const double* HP, double** ERP)
{
    const double tfw  = TF[3];
    const double rmi0 = RMI[0];
    const double rmi1 = RMI[1];
    const double hc = FPAI*Gammainv*sqrt(Gammainv)*cc;

    for (int i = 0; i < np; ++i) {
        const double Roo = 1.0 / (Zpi[i] + Gammainv);
        const double pcx = px[i] - cx;
        const double pcy = py[i] - cy;
        const double pcz = pz[i] - cz;
        const double t = Roo * (pcx*pcx + pcy*pcy + pcz*pcz);

        double f1, f2, f3;
        if (t <= tfw) {
            const int it  = (int)((t+0.015) * D33);
            const double dt  = 0.03*it - t;
            f3  = ((((ADAT[23][it]  *dt + ADAT[22][it]) * dt
                     + ADAT[21][it])*dt + ADAT[20][it]) * dt
                   + ADAT[19][it])*dt + ADAT[18][it];
            const double eed = ((((GA[5]         *dt + GA[4]) * dt
                                  + GA[3])*dt + GA[2]) * dt
                                +   GA[1])*dt + GA[0];
            const double ee  = EDAT[it] * eed;
            const double t2  = 2.0*t;
            f2  = (t2*f3 + ee)*rmi1;
            f1  = (t2*f2 + ee)*rmi0;
        } else {
            const double tinv = 1.0/t;
            f1   = 0.25*tinv*sqrt(M_PI * tinv);
            f2   = 1.5*f1*tinv;
            f3   = 2.5*f2*tinv;
        }

        const double hpq   = hc*sqrt(Roo)*HP[i];
        const double sss1  = hpq*f1;
        const double sss2  = hpq*f2;
        const double sss3  = hpq*f3;
        // [ PS | S ]
        const double Rooz  = Zpi[i]*Roo;
        const double pgx   = Rooz*pcx;
        const double pgy   = Rooz*pcy;
        const double pgz   = Rooz*pcz;
        const double xss1  = pax[i]*sss1 - pgx*sss2;
        const double yss1  = pay[i]*sss1 - pgy*sss2;
        const double zss1  = paz[i]*sss1 - pgz*sss2;
        const double xss2  = pax[i]*sss2 - pgx*sss3;
        const double yss2  = pay[i]*sss2 - pgy*sss3;
        const double zss2  = paz[i]*sss2 - pgz*sss3;
        // [ SP | S ]
        const double sxs1  = pbx[i]*sss1 - pgx*sss2;
        const double sys1  = pby[i]*sss1 - pgy*sss2;
        const double szs1  = pbz[i]*sss1 - pgz*sss2;
        // [ PP | S ]
        const double Zp5   = 0.5*Zpi[i];
        const double Zps12 = Zp5*(sss1 - Rooz*sss2);
        const double xxs1  = pbx[i]*xss1 - pgx*xss2 + Zps12;
        const double xys1  = pby[i]*xss1 - pgy*xss2;
        const double xzs1  = pbz[i]*xss1 - pgz*xss2;
        const double yxs1  = pbx[i]*yss1 - pgx*yss2;
        const double yys1  = pby[i]*yss1 - pgy*yss2 + Zps12;
        const double yzs1  = pbz[i]*yss1 - pgz*yss2;
        const double zxs1  = pbx[i]*zss1 - pgx*zss2;
        const double zys1  = pby[i]*zss1 - pgy*zss2;
        const double zzs1  = pbz[i]*zss1 - pgz*zss2 + Zps12;
        // [ PP | P ]
        const double RooG  = Roo*Gammainv;
        const double Gcx   = RooG*pcx;
        const double Gcy   = RooG*pcy;
        const double Gcz   = RooG*pcz;
        const double RZGp5 = 0.5*Rooz*Gammainv;
        const double Rxss1 = RZGp5*xss1;
        const double Ryss1 = RZGp5*yss1;
        const double Rzss1 = RZGp5*zss1;
        const double Rsxs1 = RZGp5*sxs1;
        const double Rsys1 = RZGp5*sys1;
        const double Rszs1 = RZGp5*szs1;

        ERP[i][ 0] = Gcx*xxs1 + Rsxs1 + Rxss1;
        ERP[i][ 1] = Gcy*xxs1;
        ERP[i][ 2] = Gcz*xxs1;
        ERP[i][ 3] = Gcx*xys1 + Rsys1;
        ERP[i][ 4] = Gcy*xys1         + Rxss1;
        ERP[i][ 5] = Gcz*xys1;
        ERP[i][ 6] = Gcx*xzs1 + Rszs1;
        ERP[i][ 7] = Gcy*xzs1;
        ERP[i][ 8] = Gcz*xzs1         + Rxss1;
        ERP[i][ 9] = Gcx*yxs1         + Ryss1;
        ERP[i][10] = Gcy*yxs1 + Rsxs1;
        ERP[i][11] = Gcz*yxs1;
        ERP[i][12] = Gcx*yys1;
        ERP[i][13] = Gcy*yys1 + Rsys1 + Ryss1;
        ERP[i][14] = Gcz*yys1;
        ERP[i][15] = Gcx*yzs1;
        ERP[i][16] = Gcy*yzs1 + Rszs1;
        ERP[i][17] = Gcz*yzs1         + Ryss1;
        ERP[i][18] = Gcx*zxs1         + Rzss1;
        ERP[i][19] = Gcy*zxs1;
        ERP[i][20] = Gcz*zxs1 + Rsxs1;
        ERP[i][21] = Gcx*zys1;
        ERP[i][22] = Gcy*zys1         + Rzss1;
        ERP[i][23] = Gcz*zys1 + Rsys1;
        ERP[i][24] = Gcx*zzs1;
        ERP[i][25] = Gcy*zzs1;
        ERP[i][26] = Gcz*zzs1 + Rszs1 + Rzss1;
    }
}


/*****************************************************************************/
/*                                       */
/*  Primitive 3-center ERI evaluation [ AB | C ] for type [ DS | S ]         */
/*                                       */
/*  np      ; number of first primitive pair data            */
/*  Zpi     ; 1/(Alpha+Beta) from orbital exponents Alpha, Beta  */
/*  px, py, pz  ; (Alpha*ax * Beta*bx) / (Alpha+Beta) , ...          */
/*  pax, pay, paz   ; px-ax, py-ay, pz-az                    */
/*  HP      ; Zeta**(-3/2) * exp(-(Alpha*Beta/Zeta) * (A-B)**2   */
/*  Gammainv    ; inverse of orbital exponent Gamma for last GTO     */
/*  cx, cy, cz  ; coordinates of last GTO C              */
/*  cc      ; coefficient of last GTO C              */
/*  ERP     ; output primitive ERI                   */
/*              ERP[0,1,2,...,np-1][0] for [ DxyS | S ]      */
/*              ERP[0,1,2,...,np-1][1] for [ DxzS | S ]      */
/*              ERP[0,1,2,...,np-1][2] for [ DyzS | S ]      */
/*              ERP[0,1,2,...,np-1][3] for [ Dxx-yyS | S ]   */
/*              ERP[0,1,2,...,np-1][4] for [ D3zz-rrS | S ]  */
/*                                       */
/*  FPAI   = 2 * PAI**(5/2)                          */
/*  SQR3I  = 1/sqrt(3)                           */
/*                                       */
/*****************************************************************************/
void DfEri::eriDSS(const int np,
                   const double* px, const double* py, const double* pz,
                   const double* pax, const double* pay, const double* paz,
                   const double cx, const double cy, const double cz,
                   const double cc, const double Gammainv, const double* Zpi, const double* HP, double** ERP)
{
    const double tfw  = TF[2];
    const double rmi0 = RMI[0];
    const double hc = FPAI*Gammainv*sqrt(Gammainv)*cc;

    for (int i = 0; i < np; ++i) {
        const double Roo = 1.0 / (Zpi[i] + Gammainv);
        const double pcx = px[i] - cx;
        const double pcy = py[i] - cy;
        const double pcz = pz[i] - cz;
        const double t = Roo * (pcx*pcx + pcy*pcy + pcz*pcz);

        double f0, f1, f2;
        if (t <= tfw) {
            const int it  = (int)((t+0.015) * D33);
            const double dt  = 0.03*it - t;
            f2  = ((((ADAT[17][it]  *dt + ADAT[16][it]) * dt
                     + ADAT[15][it])*dt + ADAT[14][it]) * dt
                   + ADAT[13][it])*dt + ADAT[12][it];
            const double eed = ((((GA[5]         *dt + GA[4]) * dt
                                  + GA[3])*dt + GA[2]) * dt
                                +   GA[1])*dt + GA[0];
            const double ee  = EDAT[it] * eed;
            const double t2  = 2.0*t;
            f1  = (t2*f2 + ee)*rmi0;
            f0  = t2*f1 + ee;
        } else {
            const double tinv = 1.0/t;
            f0   = 0.5*sqrt(M_PI * tinv);
            f1   = 0.5*f0*tinv;
            f2   = 1.5*f1*tinv;
        }

        const double hpq   = hc*sqrt(Roo)*HP[i];
        const double sss0  = hpq*f0;
        const double sss1  = hpq*f1;
        const double sss2  = hpq*f2;
        // [ PS | S ]
        const double Rooz  = Zpi[i]*Roo;
        const double Gpx   = -Rooz*pcx;
        const double Gpy   = -Rooz*pcy;
        const double Gpz   = -Rooz*pcz;
        const double xss0  = pax[i]*sss0 + Gpx*sss1;
        const double yss0  = pay[i]*sss0 + Gpy*sss1;
        const double zss0  = paz[i]*sss0 + Gpz*sss1;
        const double xss1  = pax[i]*sss1 + Gpx*sss2;
        const double yss1  = pay[i]*sss1 + Gpy*sss2;
        const double zss1  = paz[i]*sss1 + Gpz*sss2;
        // [ DS | S ]   xy, xz, yz, xx-yy, 3zz-rr
        const double xxssw = pax[i]*xss0 + Gpx*xss1;
        const double yyssw = pay[i]*yss0 + Gpy*yss1;

        ERP[i][0] = pay[i]*xss0 + Gpy*xss1;
        ERP[i][1] = paz[i]*xss0 + Gpz*xss1;
        ERP[i][2] = paz[i]*yss0 + Gpz*yss1;
        ERP[i][3] = 0.5*(xxssw - yyssw);
        ERP[i][4] = SQR3I*(paz[i]*zss0 + Gpz*zss1
                           - 0.5*(xxssw + yyssw));
    }
}


/*****************************************************************************/
/*                                       */
/*  Primitive 3-center ERI evaluation [ AB | C ] for type [ SS | D ]         */
/*                                       */
/*  np      ; number of first primitive pair data            */
/*  Zpi     ; 1/(Alpha+Beta) from orbital exponents Alpha, Beta  */
/*  px, py, pz  ; (Alpha*ax * Beta*bx) / (Alpha+Beta) , ...          */
/*  HP      ; Zeta**(-3/2) * exp(-(Alpha*Beta/Zeta) * (A-B)**2   */
/*  Gammainv    ; inverse of orbital exponent Gamma for last GTO     */
/*  cx, cy, cz  ; coordinates of last GTO C              */
/*  cc      ; coefficient of last GTO C              */
/*  ERP     ; output primitive ERI                   */
/*              ERP[0,1,2,...,np-1][0] for [ SS | Dxy ]      */
/*              ERP[0,1,2,...,np-1][1] for [ SS | Dxz ]      */
/*              ERP[0,1,2,...,np-1][2] for [ SS | Dyz ]      */
/*              ERP[0,1,2,...,np-1][3] for [ SS | Dxx-yy ]   */
/*              ERP[0,1,2,...,np-1][4] for [ SS | D3zz-rr ]  */
/*                                       */
/*  FPAI   = 2 * PAI**(5/2)                          */
/*  SQR3I  = 1/sqrt(3)                           */
/*                                       */
/*****************************************************************************/
void DfEri::eriSSD(const int np,
                   const double* px, const double* py, const double* pz,
                   const double cx, const double cy, const double cz,
                   const double cc, const double Gammainv, const double* Zpi, const double* HP, double** ERP)
{
    const double tfw  = TF[2];
    //const double rmi0 = RMI[0];
    const double hc   = FPAI*Gammainv*sqrt(Gammainv)*cc;

    for (int i = 0; i < np; ++i) {
        const double Roo = 1.0 / (Zpi[i] + Gammainv);
        const double pcx = px[i] - cx;
        const double pcy = py[i] - cy;
        const double pcz = pz[i] - cz;
        const double t   = Roo * (pcx*pcx + pcy*pcy + pcz*pcz);

        double f2;
        if (t <= tfw) {
            const int it  = (int)((t+0.015) * D33);
            const double dt  = 0.03*it - t;
            f2  = ((((ADAT[17][it]  *dt + ADAT[16][it]) * dt
                     + ADAT[15][it])*dt + ADAT[14][it]) * dt
                   + ADAT[13][it])*dt + ADAT[12][it];
        } else {
            const double tinv = 1.0/t;
            f2   = 0.375*sqrt(M_PI * tinv)*tinv*tinv;
        }

        const double sss2   = hc*sqrt(Roo)*HP[i]*f2;
        // [ SS | P ]
        const double RooG  = Roo*Gammainv;
        const double RRf2  = RooG*RooG*sss2;
        const double Rssx1 = pcx*RRf2;
        const double Rssy1 = pcy*RRf2;
        const double Rssz1 = pcz*RRf2;
        // [ SS | D ]   xy, xz, yz, xx-yy, 3zz-rr
        const double ssxxw = pcx*Rssx1;
        const double ssyyw = pcy*Rssy1;

        ERP[i][0] = pcx*Rssy1;
        ERP[i][1] = pcx*Rssz1;
        ERP[i][2] = pcy*Rssz1;
        ERP[i][3] = 0.5*(ssxxw - ssyyw);
        ERP[i][4] = SQR3I*(pcz*Rssz1
                           - 0.5*(ssxxw + ssyyw));
    }
}


/*****************************************************************************/
/*                                       */
/*  Primitive 3-center ERI evaluation [ AB | C ] for type [ DP | S ]         */
/*                                       */
/*  np      ; number of first primitive pair data            */
/*  Zpi     ; 1/(Alpha+Beta) from orbital exponents Alpha, Beta  */
/*  px, py, pz  ; (Alpha*ax * Beta*bx) / (Alpha+Beta) , ...          */
/*  pax, pay, paz   ; px-ax, py-ay, pz-az                    */
/*  pbx, pby, pbz   ; px-bx, py-by, pz-bz                    */
/*  HP      ; Zeta**(-3/2) * exp(-(Alpha*Beta/Zeta) * (A-B)**2   */
/*  Gammainv    ; inverse of orbital exponent Gamma for last GTO     */
/*  cx, cy, cz  ; coordinates of last GTO C              */
/*  cc      ; coefficient of last GTO C              */
/*  ERP     ; output primitive ERI                   */
/*              ERP[0,1,2,...,np-1][0] for [ DxyPx | S ]     */
/*              ERP[0,1,2,...,np-1][1] for [ DxyPy | S ]     */
/*              ERP[0,1,2,...,np-1][2] for [ DxyPz | S ]     */
/*              ERP[0,1,2,...,np-1][3] for [ DxzPx | S ]     */
/*              ERP[0,1,2,...,np-1][4] for [ DxzPy | S ]     */
/*              ERP[0,1,2,...,np-1][5] for [ DxzPz | S ]     */
/*              ERP[0,1,2,...,np-1][6] for [ DyzPx | S ]     */
/*              ERP[0,1,2,...,np-1][7] for [ DyzPy | S ]     */
/*              ERP[0,1,2,...,np-1][8] for [ DyzPz | S ]     */
/*              ERP[0,1,2,...,np-1][9] for [ Dxx-yyPx | S ]  */
/*                  .............                */
/*              ERP[0,1,2,...,np-1][12] for [ D3zz-rrPx | S ]*/
/*                  .............                */
/*                                       */
/*  FPAI   = 2 * PAI**(5/2)                          */
/*  SQR3I  = 1/sqrt(3)                           */
/*  SQR3I2 = 2/sqrt(3)                           */
/*                                       */
/*****************************************************************************/
void DfEri::eriDPS(const int np,
                   const double* px, const double* py, const double* pz,
                   const double* pax, const double* pay, const double* paz,
                   const double* pbx, const double* pby, const double* pbz,
                   const double cx, const double cy, const double cz,
                   const double cc, const double Gammainv, const double* Zpi, const double* HP, double** ERP)
{
    const double tfw  = TF[3];
    const double rmi0 = RMI[0];
    const double rmi1 = RMI[1];
    const double hc   = FPAI*Gammainv*sqrt(Gammainv)*cc;

    for (int i = 0; i < np; ++i) {
        const double Roo = 1.0 / (Zpi[i] + Gammainv);
        const double pcx = px[i] - cx;
        const double pcy = py[i] - cy;
        const double pcz = pz[i] - cz;
        const double t   = Roo * (pcx*pcx + pcy*pcy + pcz*pcz);

        double f0, f1, f2, f3;
        if (t <= tfw) {
            const int it  = (int)((t+0.015) * D33);
            const double dt  = 0.03*it - t;
            f3  = ((((ADAT[23][it]  *dt + ADAT[22][it]) * dt
                     + ADAT[21][it])*dt + ADAT[20][it]) * dt
                   + ADAT[19][it])*dt + ADAT[18][it];
            const double eed = ((((GA[5]         *dt + GA[4]) * dt
                                  + GA[3])*dt + GA[2]) * dt
                                +   GA[1])*dt + GA[0];
            const double ee  = EDAT[it] * eed;
            const double t2  = 2.0*t;
            f2  = (t2*f3 + ee)*rmi1;
            f1  = (t2*f2 + ee)*rmi0;
            f0  = t2*f1 + ee;
        } else {
            const double tinv = 1.0/t;
            f0   = 0.5*sqrt(M_PI * tinv);
            f1   = 0.5*f0*tinv;
            f2   = 1.5*f1*tinv;
            f3   = 2.5*f2*tinv;
        }

        const double hpq   = hc*sqrt(Roo)*HP[i];
        const double sss0  = hpq*f0;
        const double sss1  = hpq*f1;
        const double sss2  = hpq*f2;
        const double sss3  = hpq*f3;
        // [ PS | S ]
        const double Rooz  = Zpi[i]*Roo;
        const double Gpx   = -Rooz*pcx;
        const double Gpy   = -Rooz*pcy;
        const double Gpz   = -Rooz*pcz;
        const double xss0  = pax[i]*sss0 + Gpx*sss1;
        const double yss0  = pay[i]*sss0 + Gpy*sss1;
        const double zss0  = paz[i]*sss0 + Gpz*sss1;
        const double xss1  = pax[i]*sss1 + Gpx*sss2;
        const double yss1  = pay[i]*sss1 + Gpy*sss2;
        const double zss1  = paz[i]*sss1 + Gpz*sss2;
        const double xss2  = pax[i]*sss2 + Gpx*sss3;
        const double yss2  = pay[i]*sss2 + Gpy*sss3;
        const double zss2  = paz[i]*sss2 + Gpz*sss3;
        // [ DS | S ]   xy, xz, yz, xx-yy, 3zz-rr
        const double xxssw0 = pax[i]*xss0 + Gpx*xss1;
        const double yyssw0 = pay[i]*yss0 + Gpy*yss1;
        const double xxssw1 = pax[i]*xss1 + Gpx*xss2;
        const double yyssw1 = pay[i]*yss1 + Gpy*yss2;
        const double ass0   = pay[i]*xss0 + Gpy*xss1;
        const double bss0   = paz[i]*xss0 + Gpz*xss1;
        const double css0   = paz[i]*yss0 + Gpz*yss1;
        const double dss0   = 0.5*(xxssw0 - yyssw0);
        const double ess0   = SQR3I*(paz[i]*zss0 + Gpz*zss1
                                     - 0.5*(xxssw0 + yyssw0));
        const double ass1   = pay[i]*xss1 + Gpy*xss2;
        const double bss1   = paz[i]*xss1 + Gpz*xss2;
        const double css1   = paz[i]*yss1 + Gpz*yss2;
        const double dss1   = 0.5*(xxssw1 - yyssw1);
        const double ess1   = SQR3I*(paz[i]*zss1 + Gpz*zss2
                                     - 0.5*(xxssw1 + yyssw1));
        // [ DP | S ]
        const double Zp5   = 0.5*Zpi[i];
        const double Zxs01 = Zp5*(xss0 - Rooz*xss1);
        const double Zys01 = Zp5*(yss0 - Rooz*yss1);
        const double Zzs01 = Zp5*(zss0 - Rooz*zss1);

        ERP[i][ 0] = pbx[i]*ass0 + Gpx*ass1 + Zys01;
        ERP[i][ 1] = pby[i]*ass0 + Gpy*ass1 + Zxs01;
        ERP[i][ 2] = pbz[i]*ass0 + Gpz*ass1;
        ERP[i][ 3] = pbx[i]*bss0 + Gpx*bss1 + Zzs01;
        ERP[i][ 4] = pby[i]*bss0 + Gpy*bss1;
        ERP[i][ 5] = pbz[i]*bss0 + Gpz*bss1 + Zxs01;
        ERP[i][ 6] = pbx[i]*css0 + Gpx*css1;
        ERP[i][ 7] = pby[i]*css0 + Gpy*css1 + Zzs01;
        ERP[i][ 8] = pbz[i]*css0 + Gpz*css1 + Zys01;
        ERP[i][ 9] = pbx[i]*dss0 + Gpx*dss1 + Zxs01;
        ERP[i][10] = pby[i]*dss0 + Gpy*dss1 - Zys01;
        ERP[i][11] = pbz[i]*dss0 + Gpz*dss1;
        ERP[i][12] = pbx[i]*ess0 + Gpx*ess1 - SQR3I *Zxs01;
        ERP[i][13] = pby[i]*ess0 + Gpy*ess1 - SQR3I *Zys01;
        ERP[i][14] = pbz[i]*ess0 + Gpz*ess1 + SQR3I2*Zzs01;
    }
}


/*****************************************************************************/
/*                                       */
/*  Primitive 3-center ERI evaluation [ AB | C ] for type [ DS | P ]         */
/*                                       */
/*  np      ; number of first primitive pair data            */
/*  Zpi     ; 1/(Alpha+Beta) from orbital exponents Alpha, Beta  */
/*  px, py, pz  ; (Alpha*ax * Beta*bx) / (Alpha+Beta) , ...          */
/*  pax, pay, paz   ; px-ax, py-ay, pz-az                    */
/*  HP      ; Zeta**(-3/2) * exp(-(Alpha*Beta/Zeta) * (A-B)**2   */
/*  Gammainv    ; inverse of orbital exponent Gamma for last GTO     */
/*  cx, cy, cz  ; coordinates of last GTO C              */
/*  cc      ; coefficient of last GTO C              */
/*  ERP     ; output primitive ERI                   */
/*              ERP[0,1,2,...,np-1][0] for [ DxyS | Px ]     */
/*              ERP[0,1,2,...,np-1][1] for [ DxyS | Py ]     */
/*              ERP[0,1,2,...,np-1][2] for [ DxyS | Pz ]     */
/*              ERP[0,1,2,...,np-1][3] for [ DxzS | Px ]     */
/*              ERP[0,1,2,...,np-1][4] for [ DxzS | Py ]     */
/*              ERP[0,1,2,...,np-1][5] for [ DxzS | Pz ]     */
/*              ERP[0,1,2,...,np-1][6] for [ DyzS | Px ]     */
/*              ERP[0,1,2,...,np-1][7] for [ DyzS | Py ]     */
/*              ERP[0,1,2,...,np-1][8] for [ DyzS | Pz ]     */
/*              ERP[0,1,2,...,np-1][9] for [ Dxx-yyS | Px ]  */
/*                  .............                */
/*              ERP[0,1,2,...,np-1][12] for [ D3zz-rrS | Px ]*/
/*                  .............                */
/*                                       */
/*  FPAI   = 2 * PAI**(5/2)                          */
/*  SQR3I  = 1/sqrt(3)                           */
/*  SQR3I2 = 2/sqrt(3)                           */
/*                                       */
/*****************************************************************************/
void DfEri::eriDSP(const int np,
                   const double* px, const double* py, const double* pz,
                   const double* pax, const double* pay, const double* paz,
                   const double cx, const double cy, const double cz,
                   const double cc, const double Gammainv, const double* Zpi, const double* HP, double** ERP)
{
    const double tfw  = TF[3];
    const double rmi0 = RMI[0];
    const double rmi1 = RMI[1];
    const double hc   = FPAI*Gammainv*sqrt(Gammainv)*cc;

    for (int i = 0; i < np; ++i) {
        const double Roo = 1.0 / (Zpi[i] + Gammainv);
        const double pcx = px[i] - cx;
        const double pcy = py[i] - cy;
        const double pcz = pz[i] - cz;
        const double t   = Roo * (pcx*pcx + pcy*pcy + pcz*pcz);

        double f0, f1, f2, f3;
        if (t <= tfw) {
            const int it  = (int)((t+0.015) * D33);
            const double dt  = 0.03*it - t;
            f3  = ((((ADAT[23][it]  *dt + ADAT[22][it]) * dt
                     + ADAT[21][it])*dt + ADAT[20][it]) * dt
                   + ADAT[19][it])*dt + ADAT[18][it];
            const double eed = ((((GA[5]         *dt + GA[4]) * dt
                                  + GA[3])*dt + GA[2]) * dt
                                + GA[1])*dt + GA[0];
            const double ee  = EDAT[it] * eed;
            const double t2  = 2.0*t;
            f2  = (t2*f3 + ee)*rmi1;
            f1  = (t2*f2 + ee)*rmi0;
            f0  = t2*f1 + ee;
        } else {
            const double tinv = 1.0/t;
            f0   = 0.5*sqrt(M_PI * tinv);
            f1   = 0.5*f0*tinv;
            f2   = 1.5*f1*tinv;
            f3   = 2.5*f2*tinv;
        }

        const double hpq   = hc*sqrt(Roo)*HP[i];
        //const double sss0  = hpq*f0;
        const double sss1  = hpq*f1;
        const double sss2  = hpq*f2;
        const double sss3  = hpq*f3;
        // [ PS | S ]
        const double Rooz  = Zpi[i]*Roo;
        const double RooG  = Roo*Gammainv;
        const double Gpx   = -Rooz*pcx;
        const double Gpy   = -Rooz*pcy;
        const double Gpz   = -Rooz*pcz;
        const double Gcx   =  RooG*pcx;
        const double Gcy   =  RooG*pcy;
        const double Gcz   =  RooG*pcz;
        //const double xss0  = pax[i]*sss0 + Gpx*sss1;
        //const double yss0  = pay[i]*sss0 + Gpy*sss1;
        //const double zss0  = paz[i]*sss0 + Gpz*sss1;
        const double xss1  = pax[i]*sss1 + Gpx*sss2;
        const double yss1  = pay[i]*sss1 + Gpy*sss2;
        const double zss1  = paz[i]*sss1 + Gpz*sss2;
        const double xss2  = pax[i]*sss2 + Gpx*sss3;
        const double yss2  = pay[i]*sss2 + Gpy*sss3;
        const double zss2  = paz[i]*sss2 + Gpz*sss3;
        // [ DS | S ]   xy, xz, yz, xx-yy, 3zz-rr
        //const double xxssw0 = pax[i]*xss0 + Gpx*xss1;
        //const double yyssw0 = pay[i]*yss0 + Gpy*yss1;
        const double xxssw1 = pax[i]*xss1 + Gpx*xss2;
        const double yyssw1 = pay[i]*yss1 + Gpy*yss2;
        //const double ass0   = pay[i]*xss0 + Gpy*xss1;
        //const double bss0   = paz[i]*xss0 + Gpz*xss1;
        //const double css0   = paz[i]*yss0 + Gpz*yss1;
        //const double dss0   = 0.5*(xxssw0 - yyssw0);
        //const double ess0   = SQR3I*(paz[i]*zss0 + Gpz*zss1
        //                             - 0.5*(xxssw0 + yyssw0));
        const double ass1   = pay[i]*xss1 + Gpy*xss2;
        const double bss1   = paz[i]*xss1 + Gpz*xss2;
        const double css1   = paz[i]*yss1 + Gpz*yss2;
        const double dss1   = 0.5*(xxssw1 - yyssw1);
        const double ess1   = SQR3I*(paz[i]*zss1 + Gpz*zss2
                                     - 0.5*(xxssw1 + yyssw1));
        // [ DS | P ]
        const double RzGp5 = Rooz*Gammainv*0.5;
        const double Rxss1 = RzGp5*xss1;
        const double Ryss1 = RzGp5*yss1;
        const double Rzss1 = RzGp5*zss1;

        ERP[i][ 0] = Gcx*ass1 + Ryss1;
        ERP[i][ 1] = Gcy*ass1 + Rxss1;
        ERP[i][ 2] = Gcz*ass1;
        ERP[i][ 3] = Gcx*bss1 + Rzss1;
        ERP[i][ 4] = Gcy*bss1;
        ERP[i][ 5] = Gcz*bss1 + Rxss1;
        ERP[i][ 6] = Gcx*css1;
        ERP[i][ 7] = Gcy*css1 + Rzss1;
        ERP[i][ 8] = Gcz*css1 + Ryss1;
        ERP[i][ 9] = Gcx*dss1 + Rxss1;
        ERP[i][10] = Gcy*dss1 - Ryss1;
        ERP[i][11] = Gcz*dss1;
        ERP[i][12] = Gcx*ess1 - SQR3I *Rxss1;
        ERP[i][13] = Gcy*ess1 - SQR3I *Ryss1;
        ERP[i][14] = Gcz*ess1 + SQR3I2*Rzss1;
    }
}


/*****************************************************************************/
/*                                       */
/*  Primitive 3-center ERI evaluation [ AB | C ] for type [ PS | D ]         */
/*                                       */
/*  np      ; number of first primitive pair data            */
/*  Zpi     ; 1/(Alpha+Beta) from orbital exponents Alpha, Beta  */
/*  px, py, pz  ; (Alpha*ax * Beta*bx) / (Alpha+Beta) , ...          */
/*  pax, pay, paz   ; px-ax, py-ay, pz-az                    */
/*  HP      ; Zeta**(-3/2) * exp(-(Alpha*Beta/Zeta) * (A-B)**2   */
/*  Gammainv    ; inverse of orbital exponent Gamma for last GTO     */
/*  cx, cy, cz  ; coordinates of last GTO C              */
/*  cc      ; coefficient of last GTO C              */
/*  ERP     ; output primitive ERI                   */
/*              ERP[0,1,2,...,np-1][0] for [ PxS | Dxy ]     */
/*              ERP[0,1,2,...,np-1][1] for [ PxS | Dxz ]     */
/*              ERP[0,1,2,...,np-1][2] for [ PxS | Dyz ]     */
/*              ERP[0,1,2,...,np-1][3] for [ PxS | Dxx-yy ]  */
/*              ERP[0,1,2,...,np-1][4] for [ PxS | D3zz-rr ] */
/*              ERP[0,1,2,...,np-1][5] for [ PyS | Dxy ]     */
/*              ERP[0,1,2,...,np-1][6] for [ PyS | Dxz ]     */
/*              ERP[0,1,2,...,np-1][7] for [ PyS | Dyz ]     */
/*              ERP[0,1,2,...,np-1][8] for [ PyS | Dxx-yy]   */
/*              ERP[0,1,2,...,np-1][9] for [ PyS | D3zz-rr ] */
/*                  .............                */
/*                                       */
/*  FPAI   = 2 * PAI**(5/2)                          */
/*  SQR3I  = 1/sqrt(3)                           */
/*  SQR3I2 = 2/sqrt(3)                           */
/*                                       */
/*****************************************************************************/
void DfEri::eriPSD(const int np,
                   const double* px, const double* py, const double* pz,
                   const double* pax, const double* pay, const double* paz,
                   const double cx, const double cy, const double cz,
                   const double cc, const double Gammainv, const double* Zpi, const double* HP, double** ERP)
{
    const double tfw  = TF[3];
    const double rmi0 = RMI[0];
    const double rmi1 = RMI[1];
    const double hc   = FPAI*Gammainv*sqrt(Gammainv)*cc;

    for (int i = 0; i < np; ++i) {
        const double Roo = 1.0 / (Zpi[i] + Gammainv);
        const double pcx = px[i] - cx;
        const double pcy = py[i] - cy;
        const double pcz = pz[i] - cz;
        const double t   = Roo * (pcx*pcx + pcy*pcy + pcz*pcz);

        double f0, f1, f2, f3;
        if (t <= tfw) {
            const int it  = (int)((t+0.015) * D33);
            const double dt  = 0.03*it - t;
            f3  = ((((ADAT[23][it]  *dt + ADAT[22][it]) * dt
                     + ADAT[21][it])*dt + ADAT[20][it]) * dt
                   + ADAT[19][it])*dt + ADAT[18][it];
            const double eed = ((((GA[5]         *dt + GA[4]) * dt
                                  +GA[3])*dt + GA[2]) * dt
                                +  GA[1])*dt + GA[0];
            const double ee  = EDAT[it] * eed;
            const double t2  = 2.0*t;
            f2  = (t2*f3 + ee)*rmi1;
            f1  = (t2*f2 + ee)*rmi0;
            f0  = t2*f1 + ee;
        } else {
            const double tinv = 1.0/t;
            f0   = 0.5*sqrt(M_PI * tinv);
            f1   = 0.5*f0*tinv;
            f2   = 1.5*f1*tinv;
            f3   = 2.5*f2*tinv;
        }

        const double hpq   = hc*sqrt(Roo)*HP[i];
        const double sss2  = hpq*f2;
        const double sss3  = hpq*f3;
        // [ SS | P ]
        const double RooG  = Roo*Gammainv;
        const double Gcx   = RooG*pcx;
        const double Gcy   = RooG*pcy;
        const double Gcz   = RooG*pcz;
        const double ssx1  = Gcx*sss2;
        const double ssy1  = Gcy*sss2;
        const double ssz1  = Gcz*sss2;
        const double ssx2  = Gcx*sss3;
        const double ssy2  = Gcy*sss3;
        const double ssz2  = Gcz*sss3;
        // [ SS | D ]   xy, xz, yz, xx-yy, 3zz-rr
        const double ssxxw0 = Gcx*ssx1;
        const double ssyyw0 = Gcy*ssy1;
        const double ssxxw1 = Gcx*ssx2;
        const double ssyyw1 = Gcy*ssy2;
        const double ssa0   = Gcx*ssy1;
        const double ssb0   = Gcx*ssz1;
        const double ssc0   = Gcy*ssz1;
        const double ssd0   = 0.5*(ssxxw0 - ssyyw0);
        const double sse0   = SQR3I*(Gcz*ssz1
                                     - 0.5*(ssxxw0 + ssyyw0));
        const double ssa1   = Gcx*ssy2;
        const double ssb1   = Gcx*ssz2;
        const double ssc1   = Gcy*ssz2;
        const double ssd1   = 0.5*(ssxxw1 - ssyyw1);
        const double sse1   = SQR3I*(Gcz*ssz2
                                     - 0.5*(ssxxw1 + ssyyw1));
        // [ PS | D ]
        const double Rooz  = Roo*Zpi[i];
        const double RzGp5 = Rooz*Gammainv*0.5;
        const double Gpx   = -Rooz*pcx;
        const double Gpy   = -Rooz*pcy;
        const double Gpz   = -Rooz*pcz;
        const double Rssx1 = RzGp5*ssx1;
        const double Rssy1 = RzGp5*ssy1;
        const double Rssz1 = RzGp5*ssz1;

        ERP[i][ 0] = pax[i]*ssa0 + Gpx*ssa1 + Rssy1;
        ERP[i][ 1] = pax[i]*ssb0 + Gpx*ssb1 + Rssz1;
        ERP[i][ 2] = pax[i]*ssc0 + Gpx*ssc1;
        ERP[i][ 3] = pax[i]*ssd0 + Gpx*ssd1 + Rssx1;
        ERP[i][ 4] = pax[i]*sse0 + Gpx*sse1 - SQR3I *Rssx1;
        ERP[i][ 5] = pay[i]*ssa0 + Gpy*ssa1 + Rssx1;
        ERP[i][ 6] = pay[i]*ssb0 + Gpy*ssb1;
        ERP[i][ 7] = pay[i]*ssc0 + Gpy*ssc1 + Rssz1;
        ERP[i][ 8] = pay[i]*ssd0 + Gpy*ssd1 - Rssy1;
        ERP[i][ 9] = pay[i]*sse0 + Gpy*sse1 - SQR3I *Rssy1;
        ERP[i][10] = paz[i]*ssa0 + Gpz*ssa1;
        ERP[i][11] = paz[i]*ssb0 + Gpz*ssb1 + Rssx1;
        ERP[i][12] = paz[i]*ssc0 + Gpz*ssc1 + Rssy1;
        ERP[i][13] = paz[i]*ssd0 + Gpz*ssd1;
        ERP[i][14] = paz[i]*sse0 + Gpz*sse1 + SQR3I2*Rssz1;
    }
}


/*****************************************************************************/
/*                                       */
/*  Primitive 3-center ERI evaluation [ AB | C ] for type [ DP | P ]         */
/*                                       */
/*  np      ; number of first primitive pair data            */
/*  Zpi     ; 1/(Alpha+Beta) from orbital exponents Alpha, Beta  */
/*  px, py, pz  ; (Alpha*ax * Beta*bx) / (Alpha+Beta) , ...          */
/*  pax, pay, paz   ; px-ax, py-ay, pz-az                    */
/*  pbx, pby, pbz   ; px-bx, py-by, pz-bz                    */
/*  HP      ; Zeta**(-3/2) * exp(-(Alpha*Beta/Zeta) * (A-B)**2   */
/*  Gammainv    ; inverse of orbital exponent Gamma for last GTO     */
/*  cx, cy, cz  ; coordinates of last GTO C              */
/*  cc      ; coefficient of last GTO C              */
/*  ERP     ; output primitive ERI                   */
/*              ERP[0,1,2,...,np-1][0] for [ DxyPx | Px ]    */
/*              ERP[0,1,2,...,np-1][1] for [ DxyPx | Py ]    */
/*              ERP[0,1,2,...,np-1][2] for [ DxyPx | Pz ]    */
/*              ERP[0,1,2,...,np-1][3] for [ DxyPy | Px ]    */
/*              ERP[0,1,2,...,np-1][4] for [ DxyPy | Py ]    */
/*              ERP[0,1,2,...,np-1][5] for [ DxyPy | Pz ]    */
/*              ERP[0,1,2,...,np-1][6] for [ DxyPz | Px ]    */
/*              ERP[0,1,2,...,np-1][7] for [ DxyPz | Py ]    */
/*              ERP[0,1,2,...,np-1][8] for [ DxyPz | Pz ]    */
/*                  .............                */
/*                                       */
/*  FPAI   = 2 * PAI**(5/2)                          */
/*  SQR3I  = 1/sqrt(3)                           */
/*  SQR3I2 = 2/sqrt(3)                           */
/*                                       */
/*****************************************************************************/
void DfEri::eriDPP(const int np,
                   const double* px, const double* py, const double* pz,
                   const double* pax, const double* pay, const double* paz,
                   const double* pbx, const double* pby, const double* pbz,
                   const double cx, const double cy, const double cz,
                   const double cc, const double Gammainv, const double* Zpi, const double* HP, double** ERP)
{
    const double tfw  = TF[4];
    const double rmi0 = RMI[0];
    const double rmi1 = RMI[1];
    const double rmi2 = RMI[2];
    const double hc   = FPAI*Gammainv*sqrt(Gammainv)*cc;
    const double Gp5  = 0.5*Gammainv;

    for (int i = 0; i < np; ++i) {
        const double Roo = 1.0 / (Zpi[i] + Gammainv);
        const double pcx = px[i] - cx;
        const double pcy = py[i] - cy;
        const double pcz = pz[i] - cz;
        const double t   = Roo * (pcx*pcx + pcy*pcy + pcz*pcz);

        double f0, f1, f2, f3, f4;
        if (t <= tfw) {
            const int it  = (int)((t+0.015) * D33);
            const double dt  = 0.03*it - t;
            f4  = ((((ADAT[29][it]  *dt + ADAT[28][it]) * dt
                     + ADAT[27][it])*dt + ADAT[26][it]) * dt
                   + ADAT[25][it])*dt + ADAT[24][it];
            const double eed = ((((GA[5]         *dt + GA[4]) * dt
                                  + GA[3])*dt + GA[2]) * dt
                                +   GA[1])*dt + GA[0];
            const double ee  = EDAT[it] * eed;
            const double t2  = 2.0*t;
            f3  = (t2*f4 + ee)*rmi2;
            f2  = (t2*f3 + ee)*rmi1;
            f1  = (t2*f2 + ee)*rmi0;
            f0  = t2*f1 + ee;
        } else {
            const double tinv = 1.0/t;
            f0   = 0.5*sqrt(M_PI*tinv);
            f1   = 0.5*f0*tinv;
            f2   = 1.5*f1*tinv;
            f3   = 2.5*f2*tinv;
            f4   = 3.5*f3*tinv;
        }

        const double hpq   = hc*sqrt(Roo)*HP[i];
        //const double sss0  = hpq*f0;
        const double sss1  = hpq*f1;
        const double sss2  = hpq*f2;
        const double sss3  = hpq*f3;
        const double sss4  = hpq*f4;
        // [ PS | S ]
        const double RooG  = Roo*Gammainv;
        const double Rooz  = Zpi[i]*Roo;
        const double Gpx   = -Rooz*pcx;
        const double Gpy   = -Rooz*pcy;
        const double Gpz   = -Rooz*pcz;
        const double Gcx   =  RooG*pcx;
        const double Gcy   =  RooG*pcy;
        const double Gcz   =  RooG*pcz;
        const double xss1  = pax[i]*sss1 + Gpx*sss2;
        const double yss1  = pay[i]*sss1 + Gpy*sss2;
        const double zss1  = paz[i]*sss1 + Gpz*sss2;
        const double xss2  = pax[i]*sss2 + Gpx*sss3;
        const double yss2  = pay[i]*sss2 + Gpy*sss3;
        const double zss2  = paz[i]*sss2 + Gpz*sss3;
        const double xss3  = pax[i]*sss3 + Gpx*sss4;
        const double yss3  = pay[i]*sss3 + Gpy*sss4;
        const double zss3  = paz[i]*sss3 + Gpz*sss4;
        // [ PP | S ]
        const double Zp5   = Zpi[i]*0.5;
        const double Zps12 = Zp5*(sss1 - Rooz*sss2);
        const double xxs1  = pbx[i]*xss1 + Gpx*xss2 + Zps12;
        const double xys1  = pby[i]*xss1 + Gpy*xss2;
        const double xzs1  = pbz[i]*xss1 + Gpz*xss2;
        const double yxs1  = pbx[i]*yss1 + Gpx*yss2;
        const double yys1  = pby[i]*yss1 + Gpy*yss2 + Zps12;
        const double yzs1  = pbz[i]*yss1 + Gpz*yss2;
        const double zxs1  = pbx[i]*zss1 + Gpx*zss2;
        const double zys1  = pby[i]*zss1 + Gpy*zss2;
        const double zzs1  = pbz[i]*zss1 + Gpz*zss2 + Zps12;
        // [ DS | S ]   xy, xz, yz, xx-yy, 3zz-rr
        const double xxssw1 = pax[i]*xss1 + Gpx*xss2;
        const double yyssw1 = pay[i]*yss1 + Gpy*yss2;
        const double xxssw2 = pax[i]*xss2 + Gpx*xss3;
        const double yyssw2 = pay[i]*yss2 + Gpy*yss3;
        const double ass1   = pay[i]*xss1 + Gpy*xss2;
        const double bss1   = paz[i]*xss1 + Gpz*xss2;
        const double css1   = paz[i]*yss1 + Gpz*yss2;
        const double dss1   = 0.5*(xxssw1 - yyssw1);
        const double ess1   = SQR3I*(paz[i]*zss1 + Gpz*zss2
                                     - 0.5*(xxssw1 + yyssw1));
        const double ass2   = pay[i]*xss2 + Gpy*xss3;
        const double bss2   = paz[i]*xss2 + Gpz*xss3;
        const double css2   = paz[i]*yss2 + Gpz*yss3;
        const double dss2   = 0.5*(xxssw2 - yyssw2);
        const double ess2   = SQR3I*(paz[i]*zss2 + Gpz*zss3
                                     - 0.5*(xxssw2 + yyssw2));
        // [ DP | S ]
        const double Zxs12 = Zp5*(xss1 - Rooz*xss2);
        const double Zys12 = Zp5*(yss1 - Rooz*yss2);
        const double Zzs12 = Zp5*(zss1 - Rooz*zss2);
        const double axs1  = pbx[i]*ass1 + Gpx*ass2 + Zys12;
        const double ays1  = pby[i]*ass1 + Gpy*ass2 + Zxs12;
        const double azs1  = pbz[i]*ass1 + Gpz*ass2;
        const double bxs1  = pbx[i]*bss1 + Gpx*bss2 + Zzs12;
        const double bys1  = pby[i]*bss1 + Gpy*bss2;
        const double bzs1  = pbz[i]*bss1 + Gpz*bss2 + Zxs12;
        const double cxs1  = pbx[i]*css1 + Gpx*css2;
        const double cys1  = pby[i]*css1 + Gpy*css2 + Zzs12;
        const double czs1  = pbz[i]*css1 + Gpz*css2 + Zys12;
        const double dxs1  = pbx[i]*dss1 + Gpx*dss2 + Zxs12;
        const double dys1  = pby[i]*dss1 + Gpy*dss2 - Zys12;
        const double dzs1  = pbz[i]*dss1 + Gpz*dss2;
        const double exs1  = pbx[i]*ess1 + Gpx*ess2 - SQR3I *Zxs12;
        const double eys1  = pby[i]*ess1 + Gpy*ess2 - SQR3I *Zys12;
        const double ezs1  = pbz[i]*ess1 + Gpz*ess2 + SQR3I2*Zzs12;
        // [ DP | P ]
        const double RzGp5 = Rooz*Gp5;
        const double Rass1 = RzGp5*ass1;
        const double Rbss1 = RzGp5*bss1;
        const double Rcss1 = RzGp5*css1;
        const double Rdss1 = RzGp5*dss1;
        const double Ress1 = RzGp5*ess1;
        const double Rxxs1 = RzGp5*xxs1;
        const double Rxys1 = RzGp5*xys1;
        const double Rxzs1 = RzGp5*xzs1;
        const double Ryxs1 = RzGp5*yxs1;
        const double Ryys1 = RzGp5*yys1;
        const double Ryzs1 = RzGp5*yzs1;
        const double Rzxs1 = RzGp5*zxs1;
        const double Rzys1 = RzGp5*zys1;
        const double Rzzs1 = RzGp5*zzs1;

        // a = xy
        ERP[i][ 0] = Gcx*axs1 + Ryxs1 + Rass1;
        ERP[i][ 1] = Gcy*axs1 + Rxxs1;
        ERP[i][ 2] = Gcz*axs1;
        ERP[i][ 3] = Gcx*ays1 + Ryys1;
        ERP[i][ 4] = Gcy*ays1 + Rxys1 + Rass1;
        ERP[i][ 5] = Gcz*ays1;
        ERP[i][ 6] = Gcx*azs1 + Ryzs1;
        ERP[i][ 7] = Gcy*azs1 + Rxzs1;
        ERP[i][ 8] = Gcz*azs1         + Rass1;
        // b = xz
        ERP[i][ 9] = Gcx*bxs1 + Rzxs1 + Rbss1;
        ERP[i][10] = Gcy*bxs1;
        ERP[i][11] = Gcz*bxs1 + Rxxs1;
        ERP[i][12] = Gcx*bys1 + Rzys1;
        ERP[i][13] = Gcy*bys1         + Rbss1;
        ERP[i][14] = Gcz*bys1 + Rxys1;
        ERP[i][15] = Gcx*bzs1 + Rzzs1;
        ERP[i][16] = Gcy*bzs1;
        ERP[i][17] = Gcz*bzs1 + Rxzs1 + Rbss1;
        // c = yz
        ERP[i][18] = Gcx*cxs1         + Rcss1;
        ERP[i][19] = Gcy*cxs1 + Rzxs1;
        ERP[i][20] = Gcz*cxs1 + Ryxs1;
        ERP[i][21] = Gcx*cys1;
        ERP[i][22] = Gcy*cys1 + Rzys1 + Rcss1;
        ERP[i][23] = Gcz*cys1 + Ryys1;
        ERP[i][24] = Gcx*czs1;
        ERP[i][25] = Gcy*czs1 + Rzzs1;
        ERP[i][26] = Gcz*czs1 + Ryzs1 + Rcss1;
        // d = xx-y
        ERP[i][27] = Gcx*dxs1 + Rxxs1 + Rdss1;
        ERP[i][28] = Gcy*dxs1 - Ryxs1;
        ERP[i][29] = Gcz*dxs1;
        ERP[i][30] = Gcx*dys1 + Rxys1;
        ERP[i][31] = Gcy*dys1 - Ryys1 + Rdss1;
        ERP[i][32] = Gcz*dys1;
        ERP[i][33] = Gcx*dzs1 + Rxzs1;
        ERP[i][34] = Gcy*dzs1 - Ryzs1;
        ERP[i][35] = Gcz*dzs1         + Rdss1;
        // e = 3zz-
        ERP[i][36] = Gcx*exs1 - SQR3I *Rxxs1 + Ress1;
        ERP[i][37] = Gcy*exs1 - SQR3I *Ryxs1;
        ERP[i][38] = Gcz*exs1 + SQR3I2*Rzxs1;
        ERP[i][39] = Gcx*eys1 - SQR3I *Rxys1;
        ERP[i][40] = Gcy*eys1 - SQR3I *Ryys1 + Ress1;
        ERP[i][41] = Gcz*eys1 + SQR3I2*Rzys1;
        ERP[i][42] = Gcx*ezs1 - SQR3I *Rxzs1;
        ERP[i][43] = Gcy*ezs1 - SQR3I *Ryzs1;
        ERP[i][44] = Gcz*ezs1 + SQR3I2*Rzzs1 + Ress1;
    }
}


/*****************************************************************************/
/*                                       */
/*  Primitive 3-center ERI evaluation [ AB | C ] for type [ PP | D ]         */
/*                                       */
/*  np      ; number of first primitive pair data            */
/*  Zpi     ; 1/(Alpha+Beta) from orbital exponents Alpha, Beta  */
/*  px, py, pz  ; (Alpha*ax * Beta*bx) / (Alpha+Beta) , ...          */
/*  pax, pay, paz   ; px-ax, py-ay, pz-az                    */
/*  pbx, pby, pbz   ; px-bx, py-by, pz-bz                    */
/*  HP      ; Zeta**(-3/2) * exp(-(Alpha*Beta/Zeta) * (A-B)**2   */
/*  Gammainv    ; inverse of orbital exponent Gamma for last GTO     */
/*  cx, cy, cz  ; coordinates of last GTO C              */
/*  cc      ; coefficient of last GTO C              */
/*  ERP     ; output primitive ERI                   */
/*              ERP[0,1,2,...,np-1][0] for [ PxPx | Dxy ]    */
/*              ERP[0,1,2,...,np-1][1] for [ PxPx | Dxz ]    */
/*              ERP[0,1,2,...,np-1][2] for [ PxPx | Dyz ]    */
/*              ERP[0,1,2,...,np-1][3] for [ PxPx | Dxx-yy ] */
/*              ERP[0,1,2,...,np-1][4] for [ PxPx | D3zz-rr ]*/
/*              ERP[0,1,2,...,np-1][5] for [ PxPy | Dxy ]    */
/*              ERP[0,1,2,...,np-1][6] for [ PxPy | Dxz ]    */
/*              ERP[0,1,2,...,np-1][7] for [ PxPy | Dyz ]    */
/*              ERP[0,1,2,...,np-1][8] for [ PxPy | Dxx-yy ] */
/*              ERP[0,1,2,...,np-1][8] for [ PxPy | D3zz-rr ]*/
/*                  .............                */
/*                                       */
/*  FPAI   = 2 * PAI**(5/2)                          */
/*  SQR3I  = 1/sqrt(3)                           */
/*  SQR3I2 = 2/sqrt(3)                           */
/*                                       */
/*****************************************************************************/
void DfEri::eriPPD(const int np,
                   const double* px, const double* py, const double* pz,
                   const double* pax, const double* pay, const double* paz,
                   const double* pbx, const double* pby, const double* pbz,
                   const double cx, const double cy, const double cz,
                   const double cc, const double Gammainv, const double* Zpi, const double* HP, double** ERP)
{
    const double tfw  = TF[4];
    const double rmi0 = RMI[0];
    const double rmi1 = RMI[1];
    const double rmi2 = RMI[2];
    const double hc   = FPAI*Gammainv*sqrt(Gammainv)*cc;
    const double Gp5  = 0.5*Gammainv;

    for (int i = 0; i < np; ++i) {
        const double Roo = 1.0 / (Zpi[i] + Gammainv);
        const double pcx = px[i] - cx;
        const double pcy = py[i] - cy;
        const double pcz = pz[i] - cz;
        const double t   = Roo * (pcx*pcx + pcy*pcy + pcz*pcz);

        double f0, f1, f2, f3, f4;
        if (t <= tfw) {
            const int it  = (int)((t+0.015) * D33);
            const double dt  = 0.03*it - t;
            f4  = ((((ADAT[29][it]  *dt + ADAT[28][it]) * dt
                     + ADAT[27][it])*dt + ADAT[26][it]) * dt
                   + ADAT[25][it])*dt + ADAT[24][it];
            const double eed = ((((GA[5]         *dt + GA[4]) * dt
                                  + GA[3])*dt + GA[2]) * dt
                                +   GA[1])*dt + GA[0];
            const double ee  = EDAT[it] * eed;
            const double t2  = 2.0*t;
            f3  = (t2*f4 + ee)*rmi2;
            f2  = (t2*f3 + ee)*rmi1;
            f1  = (t2*f2 + ee)*rmi0;
            f0  = t2*f1 + ee;
        } else {
            const double tinv = 1.0/t;
            f0   = 0.5*sqrt(M_PI*tinv);
            f1   = 0.5*f0*tinv;
            f2   = 1.5*f1*tinv;
            f3   = 2.5*f2*tinv;
            f4   = 3.5*f3*tinv;
        }

        const double hpq   = hc*sqrt(Roo)*HP[i];
        //const double sss0  = hpq*f0;
        //const double sss1  = hpq*f1;
        const double sss2  = hpq*f2;
        const double sss3  = hpq*f3;
        const double sss4  = hpq*f4;
        // [ SS | P ]
        const double RooG  = Roo*Gammainv;
        const double Rooz  = Roo*Zpi[i];
        const double Zp5   = 0.5*Zpi[i];
        const double RzGp5 = Rooz*Gp5;
        const double Gpx   = -Rooz*pcx;
        const double Gpy   = -Rooz*pcy;
        const double Gpz   = -Rooz*pcz;
        const double Gcx   =  RooG*pcx;
        const double Gcy   =  RooG*pcy;
        const double Gcz   =  RooG*pcz;
        const double ssx1  = Gcx*sss2;
        const double ssy1  = Gcy*sss2;
        const double ssz1  = Gcz*sss2;
        const double ssx2  = Gcx*sss3;
        const double ssy2  = Gcy*sss3;
        const double ssz2  = Gcz*sss3;
        const double ssx3  = Gcx*sss4;
        const double ssy3  = Gcy*sss4;
        const double ssz3  = Gcz*sss4;
        // [ PS | P ]
        const double RGf2  = RzGp5*sss2;
        const double xsx1  = pax[i]*ssx1 + Gpx*ssx2 + RGf2;
        const double xsy1  = pax[i]*ssy1 + Gpx*ssy2;
        const double xsz1  = pax[i]*ssz1 + Gpx*ssz2;
        const double ysx1  = pay[i]*ssx1 + Gpy*ssx2;
        const double ysy1  = pay[i]*ssy1 + Gpy*ssy2 + RGf2;
        const double ysz1  = pay[i]*ssz1 + Gpy*ssz2;
        const double zsx1  = paz[i]*ssx1 + Gpz*ssx2;
        const double zsy1  = paz[i]*ssy1 + Gpz*ssy2;
        const double zsz1  = paz[i]*ssz1 + Gpz*ssz2 + RGf2;
        // [ SS | D ]   xy, xz, yz, xx-yy, 3zz-rr
        const double ssxxw0 = Gcx*ssx1;
        const double ssyyw0 = Gcy*ssy1;
        const double ssxxw1 = Gcx*ssx2;
        const double ssyyw1 = Gcy*ssy2;
        const double ssxxw2 = Gcx*ssx3;
        const double ssyyw2 = Gcy*ssy3;
        const double ssa0   = Gcx*ssy1;
        const double ssb0   = Gcx*ssz1;
        const double ssc0   = Gcy*ssz1;
        const double ssd0   = 0.5*(ssxxw0 - ssyyw0);
        const double sse0   = SQR3I*(Gcz*ssz1
                                     - 0.5*(ssxxw0 + ssyyw0));
        const double ssa1   = Gcx*ssy2;
        const double ssb1   = Gcx*ssz2;
        const double ssc1   = Gcy*ssz2;
        const double ssd1   = 0.5*(ssxxw1 - ssyyw1);
        const double sse1   = SQR3I*(Gcz*ssz2
                                     - 0.5*(ssxxw1 + ssyyw1));
        const double ssa2   = Gcx*ssy3;
        const double ssb2   = Gcx*ssz3;
        const double ssc2   = Gcy*ssz3;
        const double ssd2   = 0.5*(ssxxw2 - ssyyw2);
        const double sse2   = SQR3I*(Gcz*ssz3
                                     - 0.5*(ssxxw2 + ssyyw2));
        // [ PS | D ]
        const double Rssx1 = RzGp5*ssx1;
        const double Rssy1 = RzGp5*ssy1;
        const double Rssz1 = RzGp5*ssz1;
        const double xsa0  = pax[i]*ssa0 + Gpx*ssa1 + Rssy1;
        const double xsb0  = pax[i]*ssb0 + Gpx*ssb1 + Rssz1;
        const double xsc0  = pax[i]*ssc0 + Gpx*ssc1;
        const double xsd0  = pax[i]*ssd0 + Gpx*ssd1 + Rssx1;
        const double xse0  = pax[i]*sse0 + Gpx*sse1 - SQR3I *Rssx1;
        const double ysa0  = pay[i]*ssa0 + Gpy*ssa1 + Rssx1;
        const double ysb0  = pay[i]*ssb0 + Gpy*ssb1;
        const double ysc0  = pay[i]*ssc0 + Gpy*ssc1 + Rssz1;
        const double ysd0  = pay[i]*ssd0 + Gpy*ssd1 - Rssy1;
        const double yse0  = pay[i]*sse0 + Gpy*sse1 - SQR3I *Rssy1;
        const double zsa0  = paz[i]*ssa0 + Gpz*ssa1;
        const double zsb0  = paz[i]*ssb0 + Gpz*ssb1 + Rssx1;
        const double zsc0  = paz[i]*ssc0 + Gpz*ssc1 + Rssy1;
        const double zsd0  = paz[i]*ssd0 + Gpz*ssd1;
        const double zse0  = paz[i]*sse0 + Gpz*sse1 + SQR3I2*Rssz1;
        const double Rssx2 = RzGp5*ssx2;
        const double Rssy2 = RzGp5*ssy2;
        const double Rssz2 = RzGp5*ssz2;
        const double xsa1  = pax[i]*ssa1 + Gpx*ssa2 + Rssy2;
        const double xsb1  = pax[i]*ssb1 + Gpx*ssb2 + Rssz2;
        const double xsc1  = pax[i]*ssc1 + Gpx*ssc2;
        const double xsd1  = pax[i]*ssd1 + Gpx*ssd2 + Rssx2;
        const double xse1  = pax[i]*sse1 + Gpx*sse2 - SQR3I *Rssx2;
        const double ysa1  = pay[i]*ssa1 + Gpy*ssa2 + Rssx2;
        const double ysb1  = pay[i]*ssb1 + Gpy*ssb2;
        const double ysc1  = pay[i]*ssc1 + Gpy*ssc2 + Rssz2;
        const double ysd1  = pay[i]*ssd1 + Gpy*ssd2 - Rssy2;
        const double yse1  = pay[i]*sse1 + Gpy*sse2 - SQR3I *Rssy2;
        const double zsa1  = paz[i]*ssa1 + Gpz*ssa2;
        const double zsb1  = paz[i]*ssb1 + Gpz*ssb2 + Rssx2;
        const double zsc1  = paz[i]*ssc1 + Gpz*ssc2 + Rssy2;
        const double zsd1  = paz[i]*ssd1 + Gpz*ssd2;
        const double zse1  = paz[i]*sse1 + Gpz*sse2 + SQR3I2*Rssz2;
        // [ PP | D ]
        const double Zssa01 = Zp5*(ssa0 - Rooz*ssa1);
        const double Zssb01 = Zp5*(ssb0 - Rooz*ssb1);
        const double Zssc01 = Zp5*(ssc0 - Rooz*ssc1);
        const double Zssd01 = Zp5*(ssd0 - Rooz*ssd1);
        const double Zsse01 = Zp5*(sse0 - Rooz*sse1);
        const double Rxsx1  = RzGp5*xsx1;
        const double Rxsy1  = RzGp5*xsy1;
        const double Rxsz1  = RzGp5*xsz1;
        const double Rysx1  = RzGp5*ysx1;
        const double Rysy1  = RzGp5*ysy1;
        const double Rysz1  = RzGp5*ysz1;
        const double Rzsx1  = RzGp5*zsx1;
        const double Rzsy1  = RzGp5*zsy1;
        const double Rzsz1  = RzGp5*zsz1;

        ERP[i][ 0] = pbx[i]*xsa0 + Gpx*xsa1 + Zssa01 + Rxsy1;
        ERP[i][ 1] = pbx[i]*xsb0 + Gpx*xsb1 + Zssb01 + Rxsz1;
        ERP[i][ 2] = pbx[i]*xsc0 + Gpx*xsc1 + Zssc01;
        ERP[i][ 3] = pbx[i]*xsd0 + Gpx*xsd1 + Zssd01 + Rxsx1;
        ERP[i][ 4] = pbx[i]*xse0 + Gpx*xse1 + Zsse01 - SQR3I *Rxsx1;
        ERP[i][ 5] = pby[i]*xsa0 + Gpy*xsa1          + Rxsx1;
        ERP[i][ 6] = pby[i]*xsb0 + Gpy*xsb1;
        ERP[i][ 7] = pby[i]*xsc0 + Gpy*xsc1          + Rxsz1;
        ERP[i][ 8] = pby[i]*xsd0 + Gpy*xsd1          - Rxsy1;
        ERP[i][ 9] = pby[i]*xse0 + Gpy*xse1          - SQR3I *Rxsy1;
        ERP[i][10] = pbz[i]*xsa0 + Gpz*xsa1;
        ERP[i][11] = pbz[i]*xsb0 + Gpz*xsb1          + Rxsx1;
        ERP[i][12] = pbz[i]*xsc0 + Gpz*xsc1          + Rxsy1;
        ERP[i][13] = pbz[i]*xsd0 + Gpz*xsd1;
        ERP[i][14] = pbz[i]*xse0 + Gpz*xse1          + SQR3I2*Rxsz1;

        ERP[i][15] = pbx[i]*ysa0 + Gpx*ysa1          + Rysy1;
        ERP[i][16] = pbx[i]*ysb0 + Gpx*ysb1          + Rysz1;
        ERP[i][17] = pbx[i]*ysc0 + Gpx*ysc1;
        ERP[i][18] = pbx[i]*ysd0 + Gpx*ysd1          + Rysx1;
        ERP[i][19] = pbx[i]*yse0 + Gpx*yse1          - SQR3I *Rysx1;
        ERP[i][20] = pby[i]*ysa0 + Gpy*ysa1 + Zssa01 + Rysx1;
        ERP[i][21] = pby[i]*ysb0 + Gpy*ysb1 + Zssb01;
        ERP[i][22] = pby[i]*ysc0 + Gpy*ysc1 + Zssc01 + Rysz1;
        ERP[i][23] = pby[i]*ysd0 + Gpy*ysd1 + Zssd01 - Rysy1;
        ERP[i][24] = pby[i]*yse0 + Gpy*yse1 + Zsse01 - SQR3I *Rysy1;
        ERP[i][25] = pbz[i]*ysa0 + Gpz*ysa1;
        ERP[i][26] = pbz[i]*ysb0 + Gpz*ysb1          + Rysx1;
        ERP[i][27] = pbz[i]*ysc0 + Gpz*ysc1          + Rysy1;
        ERP[i][28] = pbz[i]*ysd0 + Gpz*ysd1;
        ERP[i][29] = pbz[i]*yse0 + Gpz*yse1          + SQR3I2*Rysz1;

        ERP[i][30] = pbx[i]*zsa0 + Gpx*zsa1          + Rzsy1;
        ERP[i][31] = pbx[i]*zsb0 + Gpx*zsb1          + Rzsz1;
        ERP[i][32] = pbx[i]*zsc0 + Gpx*zsc1;
        ERP[i][33] = pbx[i]*zsd0 + Gpx*zsd1          + Rzsx1;
        ERP[i][34] = pbx[i]*zse0 + Gpx*zse1          - SQR3I *Rzsx1;
        ERP[i][35] = pby[i]*zsa0 + Gpy*zsa1          + Rzsx1;
        ERP[i][36] = pby[i]*zsb0 + Gpy*zsb1;
        ERP[i][37] = pby[i]*zsc0 + Gpy*zsc1          + Rzsz1;
        ERP[i][38] = pby[i]*zsd0 + Gpy*zsd1          - Rzsy1;
        ERP[i][39] = pby[i]*zse0 + Gpy*zse1          - SQR3I *Rzsy1;
        ERP[i][40] = pbz[i]*zsa0 + Gpz*zsa1 + Zssa01;
        ERP[i][41] = pbz[i]*zsb0 + Gpz*zsb1 + Zssb01 + Rzsx1;
        ERP[i][42] = pbz[i]*zsc0 + Gpz*zsc1 + Zssc01 + Rzsy1;
        ERP[i][43] = pbz[i]*zsd0 + Gpz*zsd1 + Zssd01;
        ERP[i][44] = pbz[i]*zse0 + Gpz*zse1 + Zsse01 + SQR3I2*Rzsz1;
    }
}


/*****************************************************************************/
/*                                       */
/*  Primitive 3-center ERI evaluation [ AB | C ] for type [ DD | S ]         */
/*                                       */
/*  np      ; number of first primitive pair data            */
/*  Zpi     ; 1/(Alpha+Beta) from orbital exponents Alpha, Beta  */
/*  px, py, pz  ; (Alpha*ax * Beta*bx) / (Alpha+Beta) , ...          */
/*  pax, pay, paz   ; px-ax, py-ay, pz-az                    */
/*  pbx, pby, pbz   ; px-bx, py-by, pz-bz                    */
/*  HP      ; Zeta**(-3/2) * exp(-(Alpha*Beta/Zeta) * (A-B)**2   */
/*  Gammainv    ; inverse of orbital exponent Gamma for last GTO     */
/*  cx, cy, cz  ; coordinates of last GTO C              */
/*  cc      ; coefficient of last GTO C              */
/*  ERP     ; output primitive ERI                   */
/*              ERP[0,1,2,...,np-1][0] for [ DxyDxy | S ]    */
/*              ERP[0,1,2,...,np-1][1] for [ DxyDxz | S ]    */
/*              ERP[0,1,2,...,np-1][2] for [ DxyDyz | S ]    */
/*              ERP[0,1,2,...,np-1][3] for [ DxyDxx-yy | S ] */
/*              ERP[0,1,2,...,np-1][4] for [ DxyD3zz-rr | S ]*/
/*              ERP[0,1,2,...,np-1][5] for [ DxzDxy | S ]    */
/*              ERP[0,1,2,...,np-1][6] for [ DxzDxz | S ]    */
/*              ERP[0,1,2,...,np-1][7] for [ DxzDyz | S ]    */
/*              ERP[0,1,2,...,np-1][8] for [ DxzDxx-yy | S ] */
/*              ERP[0,1,2,...,np-1][9] for [ DxzD3zz-rr | S ]*/
/*                  .............                */
/*                                       */
/*  FPAI   = 2 * PAI**(5/2)                          */
/*  SQR3I  = 1/sqrt(3)                           */
/*  SQR3I2 = 2/sqrt(3)                           */
/*                                       */
/*****************************************************************************/
void DfEri::eriDDS(const int np,
                   const double* px, const double* py, const double* pz,
                   const double* pax, const double* pay, const double* paz,
                   const double* pbx, const double* pby, const double* pbz,
                   const double cx, const double cy, const double cz,
                   const double cc, const double Gammainv, const double* Zpi, const double* HP, double** ERP)
{
    const double tfw  = TF[4];
    const double rmi0 = RMI[0];
    const double rmi1 = RMI[1];
    const double rmi2 = RMI[2];
    const double hc   = FPAI*Gammainv*sqrt(Gammainv)*cc;
    //const double Gp5  = 0.5*Gammainv;

    for (int i = 0; i < np; ++i) {
        const double Roo = 1.0 / (Zpi[i] + Gammainv);
        const double pcx = px[i] - cx;
        const double pcy = py[i] - cy;
        const double pcz = pz[i] - cz;
        const double t   = Roo * (pcx*pcx + pcy*pcy + pcz*pcz);

        double f0, f1, f2, f3, f4;
        if (t <= tfw) {
            const int it  = (int)((t+0.015) * D33);
            const double dt  = 0.03*it - t;
            f4  = ((((ADAT[29][it]  *dt + ADAT[28][it]) * dt
                     + ADAT[27][it])*dt + ADAT[26][it]) * dt
                   + ADAT[25][it])*dt + ADAT[24][it];
            const double eed = ((((GA[5]         *dt + GA[4]) * dt
                                  + GA[3])*dt + GA[2]) * dt
                                +   GA[1])*dt + GA[0];
            const double ee  = EDAT[it] * eed;
            const double t2  = 2.0*t;
            f3  = (t2*f4 + ee)*rmi2;
            f2  = (t2*f3 + ee)*rmi1;
            f1  = (t2*f2 + ee)*rmi0;
            f0  = t2*f1 + ee;
        } else {
            const double tinv = 1.0/t;
            f0   = 0.5*sqrt(M_PI*tinv);
            f1   = 0.5*f0*tinv;
            f2   = 1.5*f1*tinv;
            f3   = 2.5*f2*tinv;
            f4   = 3.5*f3*tinv;
        }

        const double hpq   = hc*sqrt(Roo)*HP[i];
        const double sss0  = hpq*f0;
        const double sss1  = hpq*f1;
        const double sss2  = hpq*f2;
        const double sss3  = hpq*f3;
        const double sss4  = hpq*f4;
        // [ PS | S ]
        //const double RooG  = Roo*Gammainv;
        const double Rooz  = Zpi[i]*Roo;
        const double Gpx   = -Rooz*pcx;
        const double Gpy   = -Rooz*pcy;
        const double Gpz   = -Rooz*pcz;
        //const double Gcx   =  RooG*pcx;
        //const double Gcy   =  RooG*pcy;
        //const double Gcz   =  RooG*pcz;
        const double xss0  = pax[i]*sss0 + Gpx*sss1;
        const double yss0  = pay[i]*sss0 + Gpy*sss1;
        const double zss0  = paz[i]*sss0 + Gpz*sss1;
        const double xss1  = pax[i]*sss1 + Gpx*sss2;
        const double yss1  = pay[i]*sss1 + Gpy*sss2;
        const double zss1  = paz[i]*sss1 + Gpz*sss2;
        const double xss2  = pax[i]*sss2 + Gpx*sss3;
        const double yss2  = pay[i]*sss2 + Gpy*sss3;
        const double zss2  = paz[i]*sss2 + Gpz*sss3;
        const double xss3  = pax[i]*sss3 + Gpx*sss4;
        const double yss3  = pay[i]*sss3 + Gpy*sss4;
        const double zss3  = paz[i]*sss3 + Gpz*sss4;
        // [ PP | S ]
        const double Zp5   = Zpi[i]*0.5;
        //const double RzGp5 = Rooz*Gp5;
        const double Zps01 = Zp5*(sss0 - Rooz*sss1);
        const double Zps12 = Zp5*(sss1 - Rooz*sss2);
        const double xxs0  = pbx[i]*xss0 + Gpx*xss1 + Zps01;
        const double xys0  = pby[i]*xss0 + Gpy*xss1;
        const double xzs0  = pbz[i]*xss0 + Gpz*xss1;
        const double yxs0  = pbx[i]*yss0 + Gpx*yss1;
        const double yys0  = pby[i]*yss0 + Gpy*yss1 + Zps01;
        const double yzs0  = pbz[i]*yss0 + Gpz*yss1;
        const double zxs0  = pbx[i]*zss0 + Gpx*zss1;
        const double zys0  = pby[i]*zss0 + Gpy*zss1;
        const double zzs0  = pbz[i]*zss0 + Gpz*zss1 + Zps01;
        const double xxs1  = pbx[i]*xss1 + Gpx*xss2 + Zps12;
        const double xys1  = pby[i]*xss1 + Gpy*xss2;
        const double xzs1  = pbz[i]*xss1 + Gpz*xss2;
        const double yxs1  = pbx[i]*yss1 + Gpx*yss2;
        const double yys1  = pby[i]*yss1 + Gpy*yss2 + Zps12;
        const double yzs1  = pbz[i]*yss1 + Gpz*yss2;
        const double zxs1  = pbx[i]*zss1 + Gpx*zss2;
        const double zys1  = pby[i]*zss1 + Gpy*zss2;
        const double zzs1  = pbz[i]*zss1 + Gpz*zss2 + Zps12;
        // [ DS | S ]   xy, xz, yz, xx-yy, 3zz-rr
        const double xxssw0 = pax[i]*xss0 + Gpx*xss1;
        const double yyssw0 = pay[i]*yss0 + Gpy*yss1;
        const double xxssw1 = pax[i]*xss1 + Gpx*xss2;
        const double yyssw1 = pay[i]*yss1 + Gpy*yss2;
        const double xxssw2 = pax[i]*xss2 + Gpx*xss3;
        const double yyssw2 = pay[i]*yss2 + Gpy*yss3;
        const double ass0   = pay[i]*xss0 + Gpy*xss1;
        const double bss0   = paz[i]*xss0 + Gpz*xss1;
        const double css0   = paz[i]*yss0 + Gpz*yss1;
        const double dss0   = 0.5*(xxssw0 - yyssw0);
        const double ess0   = SQR3I*(paz[i]*zss0 + Gpz*zss1
                                     - 0.5*(xxssw0 + yyssw0));
        const double ass1   = pay[i]*xss1 + Gpy*xss2;
        const double bss1   = paz[i]*xss1 + Gpz*xss2;
        const double css1   = paz[i]*yss1 + Gpz*yss2;
        const double dss1   = 0.5*(xxssw1 - yyssw1);
        const double ess1   = SQR3I*(paz[i]*zss1 + Gpz*zss2
                                     - 0.5*(xxssw1 + yyssw1));
        const double ass2   = pay[i]*xss2 + Gpy*xss3;
        const double bss2   = paz[i]*xss2 + Gpz*xss3;
        const double css2   = paz[i]*yss2 + Gpz*yss3;
        const double dss2   = 0.5*(xxssw2 - yyssw2);
        const double ess2   = SQR3I*(paz[i]*zss2 + Gpz*zss3
                                     - 0.5*(xxssw2 + yyssw2));
        // [ DP | S ]
        const double Zxs01 = Zp5*(xss0 - Rooz*xss1);
        const double Zys01 = Zp5*(yss0 - Rooz*yss1);
        const double Zzs01 = Zp5*(zss0 - Rooz*zss1);
        const double Zxs12 = Zp5*(xss1 - Rooz*xss2);
        const double Zys12 = Zp5*(yss1 - Rooz*yss2);
        const double Zzs12 = Zp5*(zss1 - Rooz*zss2);
        const double axs0  = pbx[i]*ass0 + Gpx*ass1 + Zys01;
        const double ays0  = pby[i]*ass0 + Gpy*ass1 + Zxs01;
        const double azs0  = pbz[i]*ass0 + Gpz*ass1;
        const double bxs0  = pbx[i]*bss0 + Gpx*bss1 + Zzs01;
        const double bys0  = pby[i]*bss0 + Gpy*bss1;
        const double bzs0  = pbz[i]*bss0 + Gpz*bss1 + Zxs01;
        const double cxs0  = pbx[i]*css0 + Gpx*css1;
        const double cys0  = pby[i]*css0 + Gpy*css1 + Zzs01;
        const double czs0  = pbz[i]*css0 + Gpz*css1 + Zys01;
        const double dxs0  = pbx[i]*dss0 + Gpx*dss1 + Zxs01;
        const double dys0  = pby[i]*dss0 + Gpy*dss1 - Zys01;
        const double dzs0  = pbz[i]*dss0 + Gpz*dss1;
        const double exs0  = pbx[i]*ess0 + Gpx*ess1 - SQR3I *Zxs01;
        const double eys0  = pby[i]*ess0 + Gpy*ess1 - SQR3I *Zys01;
        const double ezs0  = pbz[i]*ess0 + Gpz*ess1 + SQR3I2*Zzs01;
        const double axs1  = pbx[i]*ass1 + Gpx*ass2 + Zys12;
        const double ays1  = pby[i]*ass1 + Gpy*ass2 + Zxs12;
        const double azs1  = pbz[i]*ass1 + Gpz*ass2;
        const double bxs1  = pbx[i]*bss1 + Gpx*bss2 + Zzs12;
        const double bys1  = pby[i]*bss1 + Gpy*bss2;
        const double bzs1  = pbz[i]*bss1 + Gpz*bss2 + Zxs12;
        const double cxs1  = pbx[i]*css1 + Gpx*css2;
        const double cys1  = pby[i]*css1 + Gpy*css2 + Zzs12;
        const double czs1  = pbz[i]*css1 + Gpz*css2 + Zys12;
        const double dxs1  = pbx[i]*dss1 + Gpx*dss2 + Zxs12;
        const double dys1  = pby[i]*dss1 + Gpy*dss2 - Zys12;
        const double dzs1  = pbz[i]*dss1 + Gpz*dss2;
        const double exs1  = pbx[i]*ess1 + Gpx*ess2 - SQR3I *Zxs12;
        const double eys1  = pby[i]*ess1 + Gpy*ess2 - SQR3I *Zys12;
        const double ezs1  = pbz[i]*ess1 + Gpz*ess2 + SQR3I2*Zzs12;
        // [ DD | S ]
        const double xxs01 = Zp5*(xxs0 - Rooz*xxs1);
        const double xys01 = Zp5*(xys0 - Rooz*xys1);
        const double xzs01 = Zp5*(xzs0 - Rooz*xzs1);
        const double yxs01 = Zp5*(yxs0 - Rooz*yxs1);
        const double yys01 = Zp5*(yys0 - Rooz*yys1);
        const double yzs01 = Zp5*(yzs0 - Rooz*yzs1);
        const double zxs01 = Zp5*(zxs0 - Rooz*zxs1);
        const double zys01 = Zp5*(zys0 - Rooz*zys1);
        const double zzs01 = Zp5*(zzs0 - Rooz*zzs1);
        const double axxsw = pbx[i]*axs0 + Gpx*axs1 + yxs01;
        const double ayysw = pby[i]*ays0 + Gpy*ays1 + xys01;
        const double azzsw = pbz[i]*azs0 + Gpz*azs1;
        const double bxxsw = pbx[i]*bxs0 + Gpx*bxs1 + zxs01;
        const double byysw = pby[i]*bys0 + Gpy*bys1;
        const double bzzsw = pbz[i]*bzs0 + Gpz*bzs1 + xzs01;
        const double cxxsw = pbx[i]*cxs0 + Gpx*cxs1;
        const double cyysw = pby[i]*cys0 + Gpy*cys1 + zys01;
        const double czzsw = pbz[i]*czs0 + Gpz*czs1 + yzs01;
        const double dxxsw = pbx[i]*dxs0 + Gpx*dxs1 + xxs01;
        const double dyysw = pby[i]*dys0 + Gpy*dys1 - yys01;
        const double dzzsw = pbz[i]*dzs0 + Gpz*dzs1;
        const double exxsw = pbx[i]*exs0 + Gpx*exs1 - SQR3I *xxs01;
        const double eyysw = pby[i]*eys0 + Gpy*eys1 - SQR3I *yys01;
        const double ezzsw = pbz[i]*ezs0 + Gpz*ezs1 + SQR3I2*zzs01;

        ERP[i][ 0] = pbx[i]*ays0 + Gpx*ays1 + yys01;
        ERP[i][ 1] = pbx[i]*azs0 + Gpx*azs1 + yzs01;
        ERP[i][ 2] = pby[i]*azs0 + Gpy*azs1 + xzs01;
        ERP[i][ 3] = 0.5*(axxsw - ayysw);
        ERP[i][ 4] = SQR3I*(azzsw - 0.5*(axxsw + ayysw));
        ERP[i][ 5] = pbx[i]*bys0 + Gpx*bys1 + zys01;
        ERP[i][ 6] = pbx[i]*bzs0 + Gpx*bzs1 + zzs01;
        ERP[i][ 7] = pby[i]*bzs0 + Gpy*bzs1;
        ERP[i][ 8] = 0.5*(bxxsw - byysw);
        ERP[i][ 9] = SQR3I*(bzzsw - 0.5*(bxxsw + byysw));
        ERP[i][10] = pbx[i]*cys0 + Gpx*cys1;
        ERP[i][11] = pbx[i]*czs0 + Gpx*czs1;
        ERP[i][12] = pby[i]*czs0 + Gpy*czs1 + zzs01;
        ERP[i][13] = 0.5*(cxxsw - cyysw);
        ERP[i][14] = SQR3I*(czzsw - 0.5*(cxxsw + cyysw));
        ERP[i][15] = pbx[i]*dys0 + Gpx*dys1 + xys01;
        ERP[i][16] = pbx[i]*dzs0 + Gpx*dzs1 + xzs01;
        ERP[i][17] = pby[i]*dzs0 + Gpy*dzs1 - yzs01;
        ERP[i][18] = 0.5*(dxxsw - dyysw);
        ERP[i][19] = SQR3I*(dzzsw - 0.5*(dxxsw + dyysw));
        ERP[i][20] = pbx[i]*eys0 + Gpx*eys1 - SQR3I*xys01;
        ERP[i][21] = pbx[i]*ezs0 + Gpx*ezs1 - SQR3I*xzs01;
        ERP[i][22] = pby[i]*ezs0 + Gpy*ezs1 - SQR3I*yzs01;
        ERP[i][23] = 0.5*(exxsw - eyysw);
        ERP[i][24] = SQR3I*(ezzsw - 0.5*(exxsw + eyysw));
    }
}


/*****************************************************************************/
/*                                       */
/*  Primitive 3-center ERI evaluation [ AB | C ] for type [ DS | D ]         */
/*                                       */
/*  np      ; number of first primitive pair data            */
/*  Zpi     ; 1/(Alpha+Beta) from orbital exponents Alpha, Beta  */
/*  px, py, pz  ; (Alpha*ax * Beta*bx) / (Alpha+Beta) , ...          */
/*  pax, pay, paz   ; px-ax, py-ay, pz-az                    */
/*  HP      ; Zeta**(-3/2) * exp(-(Alpha*Beta/Zeta) * (A-B)**2   */
/*  Gammainv    ; inverse of orbital exponent Gamma for last GTO     */
/*  cx, cy, cz  ; coordinates of last GTO C              */
/*  cc      ; coefficient of last GTO C              */
/*  ERP     ; output primitive ERI                   */
/*              ERP[0,1,2,...,np-1][0] for [ DxyS | Dxy ]    */
/*              ERP[0,1,2,...,np-1][1] for [ DxyS | Dxz ]    */
/*              ERP[0,1,2,...,np-1][2] for [ DxyS | Dyz ]    */
/*              ERP[0,1,2,...,np-1][3] for [ DxyS | Dxx-yy ] */
/*              ERP[0,1,2,...,np-1][4] for [ DxyS | D3zz-rr ]*/
/*              ERP[0,1,2,...,np-1][5] for [ DxzS | Dxy ]    */
/*              ERP[0,1,2,...,np-1][6] for [ DxzS | Dxz ]    */
/*              ERP[0,1,2,...,np-1][7] for [ DxzS | Dyz ]    */
/*              ERP[0,1,2,...,np-1][8] for [ DxzS | Dxx-yy ] */
/*              ERP[0,1,2,...,np-1][9] for [ DxzS | D3zz-rr ]*/
/*                  .............                */
/*                                       */
/*  FPAI   = 2 * PAI**(5/2)                          */
/*  SQR3I  = 1/sqrt(3)                           */
/*  SQR3I2 = 2/sqrt(3)                           */
/*                                       */
/*****************************************************************************/
void DfEri::eriDSD(const int np,
                   const double* px, const double* py, const double* pz,
                   const double* pax, const double* pay, const double* paz,
                   const double cx, const double cy, const double cz,
                   const double cc, const double Gammainv, const double* Zpi, const double* HP, double** ERP)
{
    const double tfw  = TF[4];
    const double rmi0 = RMI[0];
    const double rmi1 = RMI[1];
    const double rmi2 = RMI[2];
    const double hc   = FPAI*Gammainv*sqrt(Gammainv)*cc;
    const double Gp5  = 0.5*Gammainv;

    for (int i = 0; i < np; ++i) {
        const double Roo = 1.0 / (Zpi[i] + Gammainv);
        const double pcx = px[i] - cx;
        const double pcy = py[i] - cy;
        const double pcz = pz[i] - cz;
        const double t   = Roo * (pcx*pcx + pcy*pcy + pcz*pcz);

        double f0, f1, f2, f3, f4;
        if (t <= tfw) {
            const int it  = (int)((t+0.015) * D33);
            const double dt  = 0.03*it - t;
            f4  = ((((ADAT[29][it]  *dt + ADAT[28][it]) * dt
                     + ADAT[27][it])*dt + ADAT[26][it]) * dt
                   + ADAT[25][it])*dt + ADAT[24][it];
            const double eed = ((((GA[5]         *dt + GA[4]) * dt
                                  + GA[3])*dt + GA[2]) * dt
                                + GA[1])*dt + GA[0];
            const double ee  = EDAT[it] * eed;
            const double t2  = 2.0*t;
            f3  = (t2*f4 + ee)*rmi2;
            f2  = (t2*f3 + ee)*rmi1;
            f1  = (t2*f2 + ee)*rmi0;
            f0  = t2*f1 + ee;
        } else {
            const double tinv = 1.0/t;
            f0   = 0.5*sqrt(M_PI*tinv);
            f1   = 0.5*f0*tinv;
            f2   = 1.5*f1*tinv;
            f3   = 2.5*f2*tinv;
            f4   = 3.5*f3*tinv;
        }

        const double hpq   = hc*sqrt(Roo)*HP[i];
        //const double sss0  = hpq*f0;
        //const double sss1  = hpq*f1;
        const double sss2  = hpq*f2;
        const double sss3  = hpq*f3;
        const double sss4  = hpq*f4;
        // [ PS | S ]
        const double RooG  = Roo*Gammainv;
        const double Rooz  = Zpi[i]*Roo;
        const double Gpx   = -Rooz*pcx;
        const double Gpy   = -Rooz*pcy;
        const double Gpz   = -Rooz*pcz;
        const double Gcx   =  RooG*pcx;
        const double Gcy   =  RooG*pcy;
        const double Gcz   =  RooG*pcz;
        //const double xss0  = pax[i]*sss0 + Gpx*sss1;
        //const double yss0  = pay[i]*sss0 + Gpy*sss1;
        //const double zss0  = paz[i]*sss0 + Gpz*sss1;
        //const double xss1  = pax[i]*sss1 + Gpx*sss2;
        //const double yss1  = pay[i]*sss1 + Gpy*sss2;
        //const double zss1  = paz[i]*sss1 + Gpz*sss2;
        const double xss2  = pax[i]*sss2 + Gpx*sss3;
        const double yss2  = pay[i]*sss2 + Gpy*sss3;
        const double zss2  = paz[i]*sss2 + Gpz*sss3;
        const double xss3  = pax[i]*sss3 + Gpx*sss4;
        const double yss3  = pay[i]*sss3 + Gpy*sss4;
        const double zss3  = paz[i]*sss3 + Gpz*sss4;
        // [ PS | P ]
        const double RzGp5 = Rooz*Gp5;
        const double RGf2  = RzGp5*sss2;
        const double xsx1  = Gcx*xss2 + RGf2;
        const double xsy1  = Gcy*xss2;
        const double xsz1  = Gcz*xss2;
        const double ysx1  = Gcx*yss2;
        const double ysy1  = Gcy*yss2 + RGf2;
        const double ysz1  = Gcz*yss2;
        const double zsx1  = Gcx*zss2;
        const double zsy1  = Gcy*zss2;
        const double zsz1  = Gcz*zss2 + RGf2;
        // [ DS | S ]   xy, xz, yz, xx-yy, 3zz-rr
        //const double xxssw0 = pax[i]*xss0 + Gpx*xss1;
        //const double yyssw0 = pay[i]*yss0 + Gpy*yss1;
        //const double xxssw1 = pax[i]*xss1 + Gpx*xss2;
        //const double yyssw1 = pay[i]*yss1 + Gpy*yss2;
        const double xxssw2 = pax[i]*xss2 + Gpx*xss3;
        const double yyssw2 = pay[i]*yss2 + Gpy*yss3;
        //const double ass0   = pay[i]*xss0 + Gpy*xss1;
        //const double bss0   = paz[i]*xss0 + Gpz*xss1;
        //const double css0   = paz[i]*yss0 + Gpz*yss1;
        //const double dss0   = 0.5*(xxssw0 - yyssw0);
        //const double ess0   = SQR3I*(paz[i]*zss0 + Gpz*zss1
        //                             - 0.5*(xxssw0 + yyssw0));
        //const double ass1   = pay[i]*xss1 + Gpy*xss2;
        //const double bss1   = paz[i]*xss1 + Gpz*xss2;
        //const double css1   = paz[i]*yss1 + Gpz*yss2;
        //const double dss1   = 0.5*(xxssw1 - yyssw1);
        //const double ess1   = SQR3I*(paz[i]*zss1 + Gpz*zss2
        //                             - 0.5*(xxssw1 + yyssw1));
        const double ass2   = pay[i]*xss2 + Gpy*xss3;
        const double bss2   = paz[i]*xss2 + Gpz*xss3;
        const double css2   = paz[i]*yss2 + Gpz*yss3;
        const double dss2   = 0.5*(xxssw2 - yyssw2);
        const double ess2   = SQR3I*(paz[i]*zss2 + Gpz*zss3
                                     - 0.5*(xxssw2 + yyssw2));
        // [ DS | P ]
        const double Rxss2 = RzGp5*xss2;
        const double Ryss2 = RzGp5*yss2;
        const double Rzss2 = RzGp5*zss2;
        const double asx1  = Gcx*ass2 + Ryss2;
        const double asy1  = Gcy*ass2 + Rxss2;
        const double asz1  = Gcz*ass2;
        const double bsx1  = Gcx*bss2 + Rzss2;
        const double bsy1  = Gcy*bss2;
        const double bsz1  = Gcz*bss2 + Rxss2;
        const double csx1  = Gcx*css2;
        const double csy1  = Gcy*css2 + Rzss2;
        const double csz1  = Gcz*css2 + Ryss2;
        const double dsx1  = Gcx*dss2 + Rxss2;
        const double dsy1  = Gcy*dss2 - Ryss2;
        const double dsz1  = Gcz*dss2;
        const double esx1  = Gcx*ess2 - SQR3I *Rxss2;
        const double esy1  = Gcy*ess2 - SQR3I *Ryss2;
        const double esz1  = Gcz*ess2 + SQR3I2*Rzss2;
        // [ DS | D ]
        const double Rxsx1 = RzGp5*xsx1;
        const double Rxsy1 = RzGp5*xsy1;
        const double Rxsz1 = RzGp5*xsz1;
        const double Rysx1 = RzGp5*ysx1;
        const double Rysy1 = RzGp5*ysy1;
        const double Rysz1 = RzGp5*ysz1;
        const double Rzsx1 = RzGp5*zsx1;
        const double Rzsy1 = RzGp5*zsy1;
        const double Rzsz1 = RzGp5*zsz1;
        const double asxxw = Gcx*asx1 + Rysx1;
        const double asyyw = Gcy*asy1 + Rxsy1;
        const double bsxxw = Gcx*bsx1 + Rzsx1;
        const double bsyyw = Gcy*bsy1;
        const double csxxw = Gcx*csx1;
        const double csyyw = Gcy*csy1 + Rzsy1;
        const double dsxxw = Gcx*dsx1 + Rxsx1;
        const double dsyyw = Gcy*dsy1 - Rysy1;
        const double esxxw = Gcx*esx1 - SQR3I *Rxsx1;
        const double esyyw = Gcy*esy1 - SQR3I *Rysy1;

        ERP[i][ 0] = Gcx*asy1 + Rysy1;
        ERP[i][ 1] = Gcx*asz1 + Rysz1;
        ERP[i][ 2] = Gcy*asz1 + Rxsz1;
        ERP[i][ 3] = 0.5*(asxxw - asyyw);
        ERP[i][ 4] = SQR3I*(Gcz*asz1
                            - 0.5*(asxxw + asyyw));
        ERP[i][ 5] = Gcx*bsy1 + Rzsy1;
        ERP[i][ 6] = Gcx*bsz1 + Rzsz1;
        ERP[i][ 7] = Gcy*bsz1;
        ERP[i][ 8] = 0.5*(bsxxw - bsyyw);
        ERP[i][ 9] = SQR3I*(Gcz*bsz1 + Rxsz1
                            - 0.5*(bsxxw + bsyyw));
        ERP[i][10] = Gcx*csy1;
        ERP[i][11] = Gcx*csz1;
        ERP[i][12] = Gcy*csz1 + Rzsz1;
        ERP[i][13] = 0.5*(csxxw - csyyw);
        ERP[i][14] = SQR3I*(Gcz*csz1 + Rysz1
                            - 0.5*(csxxw + csyyw));
        ERP[i][15] = Gcx*dsy1 + Rxsy1;
        ERP[i][16] = Gcx*dsz1 + Rxsz1;
        ERP[i][17] = Gcy*dsz1 - Rysz1;
        ERP[i][18] = 0.5*(dsxxw - dsyyw);
        ERP[i][19] = SQR3I*(Gcz*dsz1
                            - 0.5*(dsxxw + dsyyw));
        ERP[i][20] = Gcx*esy1 - SQR3I*Rxsy1;
        ERP[i][21] = Gcx*esz1 - SQR3I*Rxsz1;
        ERP[i][22] = Gcy*esz1 - SQR3I*Rysz1;
        ERP[i][23] = 0.5*(esxxw - esyyw);
        ERP[i][24] = SQR3I*(Gcz*esz1 + SQR3I2*Rzsz1
                            - 0.5*(esxxw + esyyw));
    }
}


/*****************************************************************************/
/*                                       */
/*  Primitive 3-center ERI evaluation [ AB | C ] for type [ DD | P ]         */
/*                                       */
/*  np      ; number of first primitive pair data            */
/*  Zpi     ; 1/(Alpha+Beta) from orbital exponents Alpha, Beta  */
/*  px, py, pz  ; (Alpha*ax * Beta*bx) / (Alpha+Beta) , ...          */
/*  pax, pay, paz   ; px-ax, py-ay, pz-az                    */
/*  pbx, pby, pbz   ; px-bx, py-by, pz-bz                    */
/*  HP      ; Zeta**(-3/2) * exp(-(Alpha*Beta/Zeta) * (A-B)**2   */
/*  Gammainv    ; inverse of orbital exponent Gamma for last GTO     */
/*  cx, cy, cz  ; coordinates of last GTO C              */
/*  cc      ; coefficient of last GTO C              */
/*  ERP     ; output primitive ERI                   */
/*              ERP[0,1,2,...,np-1][0] for [ DxyDxy | Px ]   */
/*              ERP[0,1,2,...,np-1][1] for [ DxyDxy | Py ]   */
/*              ERP[0,1,2,...,np-1][2] for [ DxyDxy | Pz ]   */
/*              ERP[0,1,2,...,np-1][3] for [ DxyDxz | Px ]   */
/*              ERP[0,1,2,...,np-1][4] for [ DxyDxz | Py ]   */
/*              ERP[0,1,2,...,np-1][5] for [ DxyDxz | Pz ]   */
/*              ERP[0,1,2,...,np-1][6] for [ DxyDyz | Px ]   */
/*              ERP[0,1,2,...,np-1][7] for [ DxyDyz | Py ]   */
/*              ERP[0,1,2,...,np-1][8] for [ DxyDyz | Pz ]   */
/*              ERP[0,1,2,...,np-1][9] for [ DxyDxx-yy | Px ]*/
/*                  .............                */
/*                                       */
/*  FPAI   = 2 * PAI**(5/2)                          */
/*  SQR3I  = 1/sqrt(3)                           */
/*  SQR3I2 = 2/sqrt(3)                           */
/*                                       */
/*****************************************************************************/
void DfEri::eriDDP(const int np,
                   const double* px, const double* py, const double* pz,
                   const double* pax, const double* pay, const double* paz,
                   const double* pbx, const double* pby, const double* pbz,
                   const double cx, const double cy, const double cz,
                   const double cc, const double Gammainv, const double* Zpi, const double* HP, double** ERP)
{
    const double tfw  = TF[5];
    const double rmi0 = RMI[0];
    const double rmi1 = RMI[1];
    const double rmi2 = RMI[2];
    const double rmi3 = RMI[3];
    const double hc   = FPAI*Gammainv*sqrt(Gammainv)*cc;
    const double Gp5  = 0.5*Gammainv;

    for (int i = 0; i < np; ++i) {
        const double Roo = 1.0 / (Zpi[i] + Gammainv);
        const double pcx = px[i] - cx;
        const double pcy = py[i] - cy;
        const double pcz = pz[i] - cz;
        const double t   = Roo * (pcx*pcx + pcy*pcy + pcz*pcz);

        double f0, f1, f2, f3, f4, f5;
        if (t <= tfw) {
            const int it  = (int)((t+0.015) * D33);
            const double dt  = 0.03*it - t;
            f5  = ((((ADAT[35][it]  *dt + ADAT[34][it]) * dt
                     + ADAT[33][it])*dt + ADAT[32][it]) * dt
                   + ADAT[31][it])*dt + ADAT[30][it];
            const double eed = ((((GA[5]         *dt + GA[4]) * dt
                                  + GA[3])*dt + GA[2]) * dt
                                + GA[1])*dt + GA[0];
            const double ee  = EDAT[it] * eed;
            const double t2  = 2.0*t;
            f4  = (t2*f5 + ee)*rmi3;
            f3  = (t2*f4 + ee)*rmi2;
            f2  = (t2*f3 + ee)*rmi1;
            f1  = (t2*f2 + ee)*rmi0;
            f0  = t2*f1 + ee;
        } else {
            const double tinv = 1.0/t;
            f0   = 0.5*sqrt(M_PI*tinv);
            f1   = 0.5*f0*tinv;
            f2   = 1.5*f1*tinv;
            f3   = 2.5*f2*tinv;
            f4   = 3.5*f3*tinv;
            f5   = 4.5*f4*tinv;
        }

        const double hpq   = hc*sqrt(Roo)*HP[i];
        //const double sss0  = hpq*f0;
        const double sss1  = hpq*f1;
        const double sss2  = hpq*f2;
        const double sss3  = hpq*f3;
        const double sss4  = hpq*f4;
        const double sss5  = hpq*f5;
        // [ PS | S ]
        const double RooG  = Roo*Gammainv;
        const double Rooz  = Zpi[i]*Roo;
        const double Zp5   = Zpi[i]*0.5;
        const double RzGp5 = Rooz*Gp5;
        const double Gpx   = -Rooz*pcx;
        const double Gpy   = -Rooz*pcy;
        const double Gpz   = -Rooz*pcz;
        const double Gcx   =  RooG*pcx;
        const double Gcy   =  RooG*pcy;
        const double Gcz   =  RooG*pcz;
        const double abx   = pbx[i] - pax[i];
        const double aby   = pby[i] - pay[i];
        const double abz   = pbz[i] - paz[i];
        //const double xss0  = pax[i]*sss0 + Gpx*sss1;
        //const double yss0  = pay[i]*sss0 + Gpy*sss1;
        //const double zss0  = paz[i]*sss0 + Gpz*sss1;
        const double xss1  = pax[i]*sss1 + Gpx*sss2;
        const double yss1  = pay[i]*sss1 + Gpy*sss2;
        const double zss1  = paz[i]*sss1 + Gpz*sss2;
        const double xss2  = pax[i]*sss2 + Gpx*sss3;
        const double yss2  = pay[i]*sss2 + Gpy*sss3;
        const double zss2  = paz[i]*sss2 + Gpz*sss3;
        const double xss3  = pax[i]*sss3 + Gpx*sss4;
        const double yss3  = pay[i]*sss3 + Gpy*sss4;
        const double zss3  = paz[i]*sss3 + Gpz*sss4;
        const double xss4  = pax[i]*sss4 + Gpx*sss5;
        const double yss4  = pay[i]*sss4 + Gpy*sss5;
        const double zss4  = paz[i]*sss4 + Gpz*sss5;
        // [ PP | S ]
        //const double Zps01 = Zp5*(sss0 - Rooz*sss1);
        const double Zps12 = Zp5*(sss1 - Rooz*sss2);
        const double Zps23 = Zp5*(sss2 - Rooz*sss3);
        //const double xxs0  = pbx[i]*xss0 + Gpx*xss1 + Zps01;
        //const double xys0  = pby[i]*xss0 + Gpy*xss1;
        //const double xzs0  = pbz[i]*xss0 + Gpz*xss1;
        //const double yxs0  = pbx[i]*yss0 + Gpx*yss1;
        //const double yys0  = pby[i]*yss0 + Gpy*yss1 + Zps01;
        //const double yzs0  = pbz[i]*yss0 + Gpz*yss1;
        //const double zxs0  = pbx[i]*zss0 + Gpx*zss1;
        //const double zys0  = pby[i]*zss0 + Gpy*zss1;
        //const double zzs0  = pbz[i]*zss0 + Gpz*zss1 + Zps01;
        const double xxs1  = pbx[i]*xss1 + Gpx*xss2 + Zps12;
        const double xys1  = pby[i]*xss1 + Gpy*xss2;
        const double xzs1  = pbz[i]*xss1 + Gpz*xss2;
        const double yxs1  = pbx[i]*yss1 + Gpx*yss2;
        const double yys1  = pby[i]*yss1 + Gpy*yss2 + Zps12;
        const double yzs1  = pbz[i]*yss1 + Gpz*yss2;
        const double zxs1  = pbx[i]*zss1 + Gpx*zss2;
        const double zys1  = pby[i]*zss1 + Gpy*zss2;
        const double zzs1  = pbz[i]*zss1 + Gpz*zss2 + Zps12;
        const double xxs2  = pbx[i]*xss2 + Gpx*xss3 + Zps23;
        const double xys2  = pby[i]*xss2 + Gpy*xss3;
        const double xzs2  = pbz[i]*xss2 + Gpz*xss3;
        const double yxs2  = pbx[i]*yss2 + Gpx*yss3;
        const double yys2  = pby[i]*yss2 + Gpy*yss3 + Zps23;
        const double yzs2  = pbz[i]*yss2 + Gpz*yss3;
        const double zxs2  = pbx[i]*zss2 + Gpx*zss3;
        const double zys2  = pby[i]*zss2 + Gpy*zss3;
        const double zzs2  = pbz[i]*zss2 + Gpz*zss3 + Zps23;
        // [ DS | S ]   xy, xz, yz, xx-yy, 3zz-rr
        //const double xxssw0 = pax[i]*xss0 + Gpx*xss1;
        //const double yyssw0 = pay[i]*yss0 + Gpy*yss1;
        const double xxssw1 = pax[i]*xss1 + Gpx*xss2;
        const double yyssw1 = pay[i]*yss1 + Gpy*yss2;
        const double xxssw2 = pax[i]*xss2 + Gpx*xss3;
        const double yyssw2 = pay[i]*yss2 + Gpy*yss3;
        const double xxssw3 = pax[i]*xss3 + Gpx*xss4;
        const double yyssw3 = pay[i]*yss3 + Gpy*yss4;
        //const double ass0   = pay[i]*xss0 + Gpy*xss1;
        //const double bss0   = paz[i]*xss0 + Gpz*xss1;
        //const double css0   = paz[i]*yss0 + Gpz*yss1;
        //const double dss0   = 0.5*(xxssw0 - yyssw0);
        //const double ess0   = SQR3I*(paz[i]*zss0 + Gpz*zss1
        //                             - 0.5*(xxssw0 + yyssw0));
        const double ass1   = pay[i]*xss1 + Gpy*xss2;
        const double bss1   = paz[i]*xss1 + Gpz*xss2;
        const double css1   = paz[i]*yss1 + Gpz*yss2;
        const double dss1   = 0.5*(xxssw1 - yyssw1);
        const double ess1   = SQR3I*(paz[i]*zss1 + Gpz*zss2
                                     - 0.5*(xxssw1 + yyssw1));
        const double ass2   = pay[i]*xss2 + Gpy*xss3;
        const double bss2   = paz[i]*xss2 + Gpz*xss3;
        const double css2   = paz[i]*yss2 + Gpz*yss3;
        const double dss2   = 0.5*(xxssw2 - yyssw2);
        const double ess2   = SQR3I*(paz[i]*zss2 + Gpz*zss3
                                     - 0.5*(xxssw2 + yyssw2));
        const double ass3   = pay[i]*xss3 + Gpy*xss4;
        const double bss3   = paz[i]*xss3 + Gpz*xss4;
        const double css3   = paz[i]*yss3 + Gpz*yss4;
        const double dss3   = 0.5*(xxssw3 - yyssw3);
        const double ess3   = SQR3I*(paz[i]*zss3 + Gpz*zss4
                                     - 0.5*(xxssw3 + yyssw3));
        // [ DS | S ]   [ DxxS | S ], [ DyyS | S ], [ DzzS | S ]
        const double x2ss1  = pax[i]*xss1 + Gpx*xss2 + Zps12;
        const double y2ss1  = pay[i]*yss1 + Gpy*yss2 + Zps12;
        const double z2ss1  = paz[i]*zss1 + Gpz*zss2 + Zps12;
        const double x2ss2  = pax[i]*xss2 + Gpx*xss3 + Zps23;
        const double y2ss2  = pay[i]*yss2 + Gpy*yss3 + Zps23;
        const double z2ss2  = paz[i]*zss2 + Gpz*zss3 + Zps23;
        // [ DP | S ]
        //const double Zxs01 = Zp5*(xss0 - Rooz*xss1);
        //const double Zys01 = Zp5*(yss0 - Rooz*yss1);
        //const double Zzs01 = Zp5*(zss0 - Rooz*zss1);
        const double Zxs12 = Zp5*(xss1 - Rooz*xss2);
        const double Zys12 = Zp5*(yss1 - Rooz*yss2);
        const double Zzs12 = Zp5*(zss1 - Rooz*zss2);
        const double Zxs23 = Zp5*(xss2 - Rooz*xss3);
        const double Zys23 = Zp5*(yss2 - Rooz*yss3);
        const double Zzs23 = Zp5*(zss2 - Rooz*zss3);
        //const double axs0  = pbx[i]*ass0 + Gpx*ass1 + Zys01;
        //const double ays0  = pby[i]*ass0 + Gpy*ass1 + Zxs01;
        //const double azs0  = pbz[i]*ass0 + Gpz*ass1;
        //const double bxs0  = pbx[i]*bss0 + Gpx*bss1 + Zzs01;
        //const double bys0  = pby[i]*bss0 + Gpy*bss1;
        //const double bzs0  = pbz[i]*bss0 + Gpz*bss1 + Zxs01;
        //const double cxs0  = pbx[i]*css0 + Gpx*css1;
        //const double cys0  = pby[i]*css0 + Gpy*css1 + Zzs01;
        //const double czs0  = pbz[i]*css0 + Gpz*css1 + Zys01;
        //const double dxs0  = pbx[i]*dss0 + Gpx*dss1 + Zxs01;
        //const double dys0  = pby[i]*dss0 + Gpy*dss1 - Zys01;
        //const double dzs0  = pbz[i]*dss0 + Gpz*dss1;
        //const double exs0  = pbx[i]*ess0 + Gpx*ess1 - SQR3I *Zxs01;
        //const double eys0  = pby[i]*ess0 + Gpy*ess1 - SQR3I *Zys01;
        //const double ezs0  = pbz[i]*ess0 + Gpz*ess1 + SQR3I2*Zzs01;
        const double axs1  = pbx[i]*ass1 + Gpx*ass2 + Zys12;
        const double ays1  = pby[i]*ass1 + Gpy*ass2 + Zxs12;
        const double azs1  = pbz[i]*ass1 + Gpz*ass2;
        const double bxs1  = pbx[i]*bss1 + Gpx*bss2 + Zzs12;
        const double bys1  = pby[i]*bss1 + Gpy*bss2;
        const double bzs1  = pbz[i]*bss1 + Gpz*bss2 + Zxs12;
        const double cxs1  = pbx[i]*css1 + Gpx*css2;
        const double cys1  = pby[i]*css1 + Gpy*css2 + Zzs12;
        const double czs1  = pbz[i]*css1 + Gpz*css2 + Zys12;
        const double dxs1  = pbx[i]*dss1 + Gpx*dss2 + Zxs12;
        const double dys1  = pby[i]*dss1 + Gpy*dss2 - Zys12;
        const double dzs1  = pbz[i]*dss1 + Gpz*dss2;
        const double exs1  = pbx[i]*ess1 + Gpx*ess2 - SQR3I *Zxs12;
        const double eys1  = pby[i]*ess1 + Gpy*ess2 - SQR3I *Zys12;
        const double ezs1  = pbz[i]*ess1 + Gpz*ess2 + SQR3I2*Zzs12;
        const double axs2  = pbx[i]*ass2 + Gpx*ass3 + Zys23;
        const double ays2  = pby[i]*ass2 + Gpy*ass3 + Zxs23;
        const double azs2  = pbz[i]*ass2 + Gpz*ass3;
        const double bxs2  = pbx[i]*bss2 + Gpx*bss3 + Zzs23;
        const double bys2  = pby[i]*bss2 + Gpy*bss3;
        const double bzs2  = pbz[i]*bss2 + Gpz*bss3 + Zxs23;
        const double cxs2  = pbx[i]*css2 + Gpx*css3;
        const double cys2  = pby[i]*css2 + Gpy*css3 + Zzs23;
        const double czs2  = pbz[i]*css2 + Gpz*css3 + Zys23;
        const double dxs2  = pbx[i]*dss2 + Gpx*dss3 + Zxs23;
        const double dys2  = pby[i]*dss2 + Gpy*dss3 - Zys23;
        const double dzs2  = pbz[i]*dss2 + Gpz*dss3;
        const double exs2  = pbx[i]*ess2 + Gpx*ess3 - SQR3I *Zxs23;
        const double eys2  = pby[i]*ess2 + Gpy*ess3 - SQR3I *Zys23;
        const double ezs2  = pbz[i]*ess2 + Gpz*ess3 + SQR3I2*Zzs23;
        // [ DP | S ]   [ DxxPx | S ], [ DyyPy | S ], [ DzzPz | S ]
        const double x2xs1  = pbx[i]*x2ss1 + Gpx*x2ss2 + 2.0*Zxs12;
        const double y2ys1  = pby[i]*y2ss1 + Gpy*y2ss2 + 2.0*Zys12;
        const double z2zs1  = pbz[i]*z2ss1 + Gpz*z2ss2 + 2.0*Zzs12;
        // [ PD | S ]
        const double xas1 = axs1 + aby*xxs1;
        const double xbs1 = bxs1 + abz*xxs1;
        const double xcs1 = azs1 + aby*xzs1;
        const double xds1 = 0.5*(x2xs1 + abx*xxs1 - ays1 - aby*xys1);
        const double xes1 = SQR3I*(bzs1 + abz*xzs1
                                   -0.5*(x2xs1 + abx*xxs1 + ays1 + aby*xys1));
        const double yas1 = ays1 + abx*yys1;
        const double ybs1 = cxs1 + abz*yxs1;
        const double ycs1 = cys1 + abz*yys1;
        const double yds1 = 0.5*(axs1 + abx*yxs1 - y2ys1 - aby*yys1);
        const double yes1 = SQR3I*(czs1 + abz*yzs1
                                   -0.5*(axs1 + abx*yxs1 + y2ys1 + aby*yys1));
        const double zas1 = bys1 + abx*zys1;
        const double zbs1 = bzs1 + abx*zzs1;
        const double zcs1 = czs1 + aby*zzs1;
        const double zds1 = 0.5*(bxs1 + abx*zxs1 - cys1 - aby*zys1);
        const double zes1 = SQR3I*(z2zs1 + abz*zzs1
                                   -0.5*(bxs1 + abx*zxs1 + cys1 + aby*zys1));
        // [ DP | P ]
        //const double Rass2 = RzGp5*ass2;
        //const double Rbss2 = RzGp5*bss2;
        //const double Rcss2 = RzGp5*css2;
        //const double Rdss2 = RzGp5*dss2;
        //const double Ress2 = RzGp5*ess2;
        //const double Rxxs2 = RzGp5*xxs2;
        //const double Rxys2 = RzGp5*xys2;
        //const double Rxzs2 = RzGp5*xzs2;
        //const double Ryxs2 = RzGp5*yxs2;
        //const double Ryys2 = RzGp5*yys2;
        //const double Ryzs2 = RzGp5*yzs2;
        //const double Rzxs2 = RzGp5*zxs2;
        //const double Rzys2 = RzGp5*zys2;
        //const double Rzzs2 = RzGp5*zzs2;
        // a = xy
        //const double axx1 = Gcx*axs2 + Ryxs2 + Rass2;
        //const double axy1 = Gcy*axs2 + Rxxs2;
        //const double axz1 = Gcz*axs2;
        //const double ayx1 = Gcx*ays2 + Ryys2;
        //const double ayy1 = Gcy*ays2 + Rxys2 + Rass2;
        //const double ayz1 = Gcz*ays2;
        //const double azx1 = Gcx*azs2 + Ryzs2;
        //const double azy1 = Gcy*azs2 + Rxzs2;
        //const double azz1 = Gcz*azs2         + Rass2;
        // b = xz
        //const double bxx1 = Gcx*bxs2 + Rzxs2 + Rbss2;
        //const double bxy1 = Gcy*bxs2;
        //const double bxz1 = Gcz*bxs2 + Rxxs2;
        //const double byx1 = Gcx*bys2 + Rzys2;
        //const double byy1 = Gcy*bys2         + Rbss2;
        //const double byz1 = Gcz*bys2 + Rxys2;
        //const double bzx1 = Gcx*bzs2 + Rzzs2;
        //const double bzy1 = Gcy*bzs2;
        //const double bzz1 = Gcz*bzs2 + Rxzs2 + Rbss2;
        // c = yz
        //const double cxx1 = Gcx*cxs2         + Rcss2;
        //const double cxy1 = Gcy*cxs2 + Rzxs2;
        //const double cxz1 = Gcz*cxs2 + Ryxs2;
        //const double cyx1 = Gcx*cys2;
        //const double cyy1 = Gcy*cys2 + Rzys2 + Rcss2;
        //const double cyz1 = Gcz*cys2 + Ryys2;
        //const double czx1 = Gcx*czs2;
        //const double czy1 = Gcy*czs2 + Rzzs2;
        //const double czz1 = Gcz*czs2 + Ryzs2 + Rcss2;
        // d = xx-yy
        //const double dxx1 = Gcx*dxs2 + Rxxs2 + Rdss2;
        //const double dxy1 = Gcy*dxs2 - Ryxs2;
        //const double dxz1 = Gcz*dxs2;
        //const double dyx1 = Gcx*dys2 + Rxys2;
        //const double dyy1 = Gcy*dys2 - Ryys2 + Rdss2;
        //const double dyz1 = Gcz*dys2;
        //const double dzx1 = Gcx*dzs2 + Rxzs2;
        //const double dzy1 = Gcy*dzs2 - Ryzs2;
        //const double dzz1 = Gcz*dzs2         + Rdss2;
        // e = 3zz-rr
        //const double exx1 = Gcx*exs2 - SQR3I *Rxxs2 + Ress2;
        //const double exy1 = Gcy*exs2 - SQR3I *Ryxs2;
        //const double exz1 = Gcz*exs2 + SQR3I2*Rzxs2;
        //const double eyx1 = Gcx*eys2 - SQR3I *Rxys2;
        //const double eyy1 = Gcy*eys2 - SQR3I *Ryys2 + Ress2;
        //const double eyz1 = Gcz*eys2 + SQR3I2*Rzys2;
        //const double ezx1 = Gcx*ezs2 - SQR3I *Rxzs2;
        //const double ezy1 = Gcy*ezs2 - SQR3I *Ryzs2;
        //const double ezz1 = Gcz*ezs2 + SQR3I2*Rzzs2 + Ress2;
        // [ DD | S ]
        const double xxs12 = Zp5*(xxs1 - Rooz*xxs2);
        const double xys12 = Zp5*(xys1 - Rooz*xys2);
        const double xzs12 = Zp5*(xzs1 - Rooz*xzs2);
        const double yxs12 = Zp5*(yxs1 - Rooz*yxs2);
        const double yys12 = Zp5*(yys1 - Rooz*yys2);
        const double yzs12 = Zp5*(yzs1 - Rooz*yzs2);
        const double zxs12 = Zp5*(zxs1 - Rooz*zxs2);
        const double zys12 = Zp5*(zys1 - Rooz*zys2);
        const double zzs12 = Zp5*(zzs1 - Rooz*zzs2);
        const double axxsw = pbx[i]*axs1 + Gpx*axs2 + yxs12;
        const double ayysw = pby[i]*ays1 + Gpy*ays2 + xys12;
        const double azzsw = pbz[i]*azs1 + Gpz*azs2;
        const double bxxsw = pbx[i]*bxs1 + Gpx*bxs2 + zxs12;
        const double byysw = pby[i]*bys1 + Gpy*bys2;
        const double bzzsw = pbz[i]*bzs1 + Gpz*bzs2 + xzs12;
        const double cxxsw = pbx[i]*cxs1 + Gpx*cxs2;
        const double cyysw = pby[i]*cys1 + Gpy*cys2 + zys12;
        const double czzsw = pbz[i]*czs1 + Gpz*czs2 + yzs12;
        const double dxxsw = pbx[i]*dxs1 + Gpx*dxs2 + xxs12;
        const double dyysw = pby[i]*dys1 + Gpy*dys2 - yys12;
        const double dzzsw = pbz[i]*dzs1 + Gpz*dzs2;
        const double exxsw = pbx[i]*exs1 + Gpx*exs2 - SQR3I *xxs12;
        const double eyysw = pby[i]*eys1 + Gpy*eys2 - SQR3I *yys12;
        const double ezzsw = pbz[i]*ezs1 + Gpz*ezs2 + SQR3I2*zzs12;

        const double aas1  = pbx[i]*ays1 + Gpx*ays2 + yys12;
        const double abs1  = pbx[i]*azs1 + Gpx*azs2 + yzs12;
        const double acs1  = pby[i]*azs1 + Gpy*azs2 + xzs12;
        const double ads1  = 0.5*(axxsw - ayysw);
        const double aes1  = SQR3I*(azzsw - 0.5*(axxsw + ayysw));
        const double bas1  = pbx[i]*bys1 + Gpx*bys2 + zys12;
        const double bbs1  = pbx[i]*bzs1 + Gpx*bzs2 + zzs12;
        const double bcs1  = pby[i]*bzs1 + Gpy*bzs2;
        const double bds1  = 0.5*(bxxsw - byysw);
        const double bes1  = SQR3I*(bzzsw - 0.5*(bxxsw + byysw));
        const double cas1  = pbx[i]*cys1 + Gpx*cys2;
        const double cbs1  = pbx[i]*czs1 + Gpx*czs2;
        const double ccs1  = pby[i]*czs1 + Gpy*czs2 + zzs12;
        const double cds1  = 0.5*(cxxsw - cyysw);
        const double ces1  = SQR3I*(czzsw - 0.5*(cxxsw + cyysw));
        const double das1  = pbx[i]*dys1 + Gpx*dys2 + xys12;
        const double dbs1  = pbx[i]*dzs1 + Gpx*dzs2 + xzs12;
        const double dcs1  = pby[i]*dzs1 + Gpy*dzs2 - yzs12;
        const double dds1  = 0.5*(dxxsw - dyysw);
        const double des1  = SQR3I*(dzzsw - 0.5*(dxxsw + dyysw));
        const double eas1  = pbx[i]*eys1 + Gpx*eys2 - SQR3I*xys12;
        const double ebs1  = pbx[i]*ezs1 + Gpx*ezs2 - SQR3I*xzs12;
        const double ecs1  = pby[i]*ezs1 + Gpy*ezs2 - SQR3I*yzs12;
        const double eds1  = 0.5*(exxsw - eyysw);
        const double ees1  = SQR3I*(ezzsw - 0.5*(exxsw + eyysw));
        // [ DD | P ]
        const double Raxs1 = RzGp5*axs1;
        const double Rays1 = RzGp5*ays1;
        const double Razs1 = RzGp5*azs1;
        const double Rbxs1 = RzGp5*bxs1;
        const double Rbys1 = RzGp5*bys1;
        const double Rbzs1 = RzGp5*bzs1;
        const double Rcxs1 = RzGp5*cxs1;
        const double Rcys1 = RzGp5*cys1;
        const double Rczs1 = RzGp5*czs1;
        const double Rdxs1 = RzGp5*dxs1;
        const double Rdys1 = RzGp5*dys1;
        const double Rdzs1 = RzGp5*dzs1;
        const double Rexs1 = RzGp5*exs1;
        const double Reys1 = RzGp5*eys1;
        const double Rezs1 = RzGp5*ezs1;
        const double Rxas1 = RzGp5*xas1;
        const double Rxbs1 = RzGp5*xbs1;
        const double Rxcs1 = RzGp5*xcs1;
        const double Rxds1 = RzGp5*xds1;
        const double Rxes1 = RzGp5*xes1;
        const double Ryas1 = RzGp5*yas1;
        const double Rybs1 = RzGp5*ybs1;
        const double Rycs1 = RzGp5*ycs1;
        const double Ryds1 = RzGp5*yds1;
        const double Ryes1 = RzGp5*yes1;
        const double Rzas1 = RzGp5*zas1;
        const double Rzbs1 = RzGp5*zbs1;
        const double Rzcs1 = RzGp5*zcs1;
        const double Rzds1 = RzGp5*zds1;
        const double Rzes1 = RzGp5*zes1;

        ERP[i][ 0] = Gcx*aas1 + Ryas1 + Rays1;
        ERP[i][ 1] = Gcy*aas1 + Rxas1 + Raxs1;
        ERP[i][ 2] = Gcz*aas1;
        ERP[i][ 3] = Gcx*abs1 + Rybs1 + Razs1;
        ERP[i][ 4] = Gcy*abs1 + Rxbs1;
        ERP[i][ 5] = Gcz*abs1         + Raxs1;
        ERP[i][ 6] = Gcx*acs1 + Rycs1;
        ERP[i][ 7] = Gcy*acs1 + Rxcs1 + Razs1;
        ERP[i][ 8] = Gcz*acs1         + Rays1;
        ERP[i][ 9] = Gcx*ads1 + Ryds1 + Raxs1;
        ERP[i][10] = Gcy*ads1 + Rxds1 - Rays1;
        ERP[i][11] = Gcz*ads1;
        ERP[i][12] = Gcx*aes1 + Ryes1 - SQR3I *Raxs1;
        ERP[i][13] = Gcy*aes1 + Rxes1 - SQR3I *Rays1;
        ERP[i][14] = Gcz*aes1         + SQR3I2*Razs1;

        ERP[i][15] = Gcx*bas1 + Rzas1 + Rbys1;
        ERP[i][16] = Gcy*bas1         + Rbxs1;
        ERP[i][17] = Gcz*bas1 + Rxas1;
        ERP[i][18] = Gcx*bbs1 + Rzbs1 + Rbzs1;
        ERP[i][19] = Gcy*bbs1;
        ERP[i][20] = Gcz*bbs1 + Rxbs1 + Rbxs1;
        ERP[i][21] = Gcx*bcs1 + Rzcs1;
        ERP[i][22] = Gcy*bcs1         + Rbzs1;
        ERP[i][23] = Gcz*bcs1 + Rxcs1 + Rbys1;
        ERP[i][24] = Gcx*bds1 + Rzds1 + Rbxs1;
        ERP[i][25] = Gcy*bds1         - Rbys1;
        ERP[i][26] = Gcz*bds1 + Rxds1;
        ERP[i][27] = Gcx*bes1 + Rzes1 - SQR3I *Rbxs1;
        ERP[i][28] = Gcy*bes1         - SQR3I *Rbys1;
        ERP[i][29] = Gcz*bes1 + Rxes1 + SQR3I2*Rbzs1;

        ERP[i][30] = Gcx*cas1         + Rcys1;
        ERP[i][31] = Gcy*cas1 + Rzas1 + Rcxs1;
        ERP[i][32] = Gcz*cas1 + Ryas1;
        ERP[i][33] = Gcx*cbs1         + Rczs1;
        ERP[i][34] = Gcy*cbs1 + Rzbs1;
        ERP[i][35] = Gcz*cbs1 + Rybs1 + Rcxs1;
        ERP[i][36] = Gcx*ccs1;
        ERP[i][37] = Gcy*ccs1 + Rzcs1 + Rczs1;
        ERP[i][38] = Gcz*ccs1 + Rycs1 + Rcys1;
        ERP[i][39] = Gcx*cds1         + Rcxs1;
        ERP[i][40] = Gcy*cds1 + Rzds1 - Rcys1;
        ERP[i][41] = Gcz*cds1 + Ryds1;
        ERP[i][42] = Gcx*ces1         - SQR3I *Rcxs1;
        ERP[i][43] = Gcy*ces1 + Rzes1 - SQR3I *Rcys1;
        ERP[i][44] = Gcz*ces1 + Ryes1 + SQR3I2*Rczs1;

        ERP[i][45] = Gcx*das1 + Rxas1 + Rdys1;
        ERP[i][46] = Gcy*das1 - Ryas1 + Rdxs1;
        ERP[i][47] = Gcz*das1;
        ERP[i][48] = Gcx*dbs1 + Rxbs1 + Rdzs1;
        ERP[i][49] = Gcy*dbs1 - Rybs1;
        ERP[i][50] = Gcz*dbs1         + Rdxs1;
        ERP[i][51] = Gcx*dcs1 + Rxcs1;
        ERP[i][52] = Gcy*dcs1 - Rycs1 + Rdzs1;
        ERP[i][53] = Gcz*dcs1         + Rdys1;
        ERP[i][54] = Gcx*dds1 + Rxds1 + Rdxs1;
        ERP[i][55] = Gcy*dds1 - Ryds1 - Rdys1;
        ERP[i][56] = Gcz*dds1;
        ERP[i][57] = Gcx*des1 + Rxes1 - SQR3I *Rdxs1;
        ERP[i][58] = Gcy*des1 - Ryes1 - SQR3I *Rdys1;
        ERP[i][59] = Gcz*des1         + SQR3I2*Rdzs1;

        ERP[i][60] = Gcx*eas1 - SQR3I *Rxas1 + Reys1;
        ERP[i][61] = Gcy*eas1 - SQR3I *Ryas1 + Rexs1;
        ERP[i][62] = Gcz*eas1 + SQR3I2*Rzas1;
        ERP[i][63] = Gcx*ebs1 - SQR3I *Rxbs1 + Rezs1;
        ERP[i][64] = Gcy*ebs1 - SQR3I *Rybs1;
        ERP[i][65] = Gcz*ebs1 + SQR3I2*Rzbs1 + Rexs1;
        ERP[i][66] = Gcx*ecs1 - SQR3I *Rxcs1;
        ERP[i][67] = Gcy*ecs1 - SQR3I *Rycs1 + Rezs1;
        ERP[i][68] = Gcz*ecs1 + SQR3I2*Rzcs1 + Reys1;
        ERP[i][69] = Gcx*eds1 - SQR3I *Rxds1 + Rexs1;
        ERP[i][70] = Gcy*eds1 - SQR3I *Ryds1 - Reys1;
        ERP[i][71] = Gcz*eds1 + SQR3I2*Rzds1;
        ERP[i][72] = Gcx*ees1 - SQR3I *(Rxes1 + Rexs1);
        ERP[i][73] = Gcy*ees1 - SQR3I *(Ryes1 + Reys1);
        ERP[i][74] = Gcz*ees1 + SQR3I2*(Rzes1 + Rezs1);
    }
}


/*****************************************************************************/
/*                                       */
/*  Primitive 3-center ERI evaluation [ AB | C ] for type [ DP | D ]         */
/*                                       */
/*  np      ; number of first primitive pair data            */
/*  Zpi     ; 1/(Alpha+Beta) from orbital exponents Alpha, Beta  */
/*  px, py, pz  ; (Alpha*ax * Beta*bx) / (Alpha+Beta) , ...          */
/*  pax, pay, paz   ; px-ax, py-ay, pz-az                    */
/*  pbx, pby, pbz   ; px-bx, py-by, pz-bz                    */
/*  HP      ; Zeta**(-3/2) * exp(-(Alpha*Beta/Zeta) * (A-B)**2   */
/*  Gammainv    ; inverse of orbital exponent Gamma for last GTO     */
/*  cx, cy, cz  ; coordinates of last GTO C              */
/*  cc      ; coefficient of last GTO C              */
/*  ERP     ; output primitive ERI                   */
/*          ERP[0,1,2,...,np-1][0] for [ DxyPx | Dxy ]       */
/*          ERP[0,1,2,...,np-1][1] for [ DxyPx | Dxz ]       */
/*          ERP[0,1,2,...,np-1][2] for [ DxyPx | Dyz ]       */
/*          ERP[0,1,2,...,np-1][3] for [ DxyPx | Dxx-yy ]        */
/*          ERP[0,1,2,...,np-1][4] for [ DxyPx | D3zz-rr ]       */
/*          ERP[0,1,2,...,np-1][5] for [ DxyPy | Dxy ]       */
/*          ERP[0,1,2,...,np-1][6] for [ DxyPy | Dxz ]       */
/*          ERP[0,1,2,...,np-1][7] for [ DxyPy | Dyz ]           */
/*          ERP[0,1,2,...,np-1][8] for [ DxyPy | Dxx-yy ]        */
/*          ERP[0,1,2,...,np-1][9] for [ DxyPy | D3zz-rr ]       */
/*                  .............                */
/*                                       */
/*  FPAI   = 2 * PAI**(5/2)                          */
/*  SQR3I  = 1/sqrt(3)                           */
/*  SQR3I2 = 2/sqrt(3)                           */
/*                                       */
/*****************************************************************************/
void DfEri::eriDPD(const int np,
                   const double* px, const double* py, const double* pz,
                   const double* pax, const double* pay, const double* paz,
                   const double* pbx, const double* pby, const double* pbz,
                   const double cx, const double cy, const double cz,
                   const double cc, const double Gammainv, const double* Zpi, const double* HP, double** ERP)
{
    const double tfw  = TF[5];
    const double rmi0 = RMI[0];
    const double rmi1 = RMI[1];
    const double rmi2 = RMI[2];
    const double rmi3 = RMI[3];
    const double hc   = FPAI*Gammainv*sqrt(Gammainv)*cc;
    const double Gp5  = 0.5*Gammainv;

    for (int i = 0; i < np; ++i) {
        const double Roo = 1.0 / (Zpi[i] + Gammainv);
        const double pcx = px[i] - cx;
        const double pcy = py[i] - cy;
        const double pcz = pz[i] - cz;
        const double t   = Roo * (pcx*pcx + pcy*pcy + pcz*pcz);

        double f0, f1, f2, f3, f4, f5;
        if (t <= tfw) {
            const int it  = (int)((t+0.015) * D33);
            const double dt  = 0.03*it - t;
            f5  = ((((ADAT[35][it]  *dt + ADAT[34][it]) * dt
                     + ADAT[33][it])*dt + ADAT[32][it]) * dt
                   + ADAT[31][it])*dt + ADAT[30][it];
            const double  eed = ((((GA[5]         *dt + GA[4]) * dt
                                   + GA[3])*dt + GA[2]) * dt
                                 + GA[1])*dt + GA[0];
            const double ee  = EDAT[it] * eed;
            const double t2  = 2.0*t;
            f4  = (t2*f5 + ee)*rmi3;
            f3  = (t2*f4 + ee)*rmi2;
            f2  = (t2*f3 + ee)*rmi1;
            f1  = (t2*f2 + ee)*rmi0;
            f0  = t2*f1 + ee;
        } else {
            const double tinv = 1.0/t;
            f0   = 0.5*sqrt(M_PI*tinv);
            f1   = 0.5*f0*tinv;
            f2   = 1.5*f1*tinv;
            f3   = 2.5*f2*tinv;
            f4   = 3.5*f3*tinv;
            f5   = 4.5*f4*tinv;
        }

        const double hpq   = hc*sqrt(Roo)*HP[i];
        //const double sss0  = hpq*f0;
        //const double sss1  = hpq*f1;
        const double sss2  = hpq*f2;
        const double sss3  = hpq*f3;
        const double sss4  = hpq*f4;
        const double sss5  = hpq*f5;
        // [ PS | S ]
        const double RooG  = Roo*Gammainv;
        const double Rooz  = Zpi[i]*Roo;
        const double Zp5   = Zpi[i]*0.5;
        const double RzGp5 = Rooz*Gp5;
        const double Gpx   = -Rooz*pcx;
        const double Gpy   = -Rooz*pcy;
        const double Gpz   = -Rooz*pcz;
        const double Gcx   =  RooG*pcx;
        const double Gcy   =  RooG*pcy;
        const double Gcz   =  RooG*pcz;
        //const double xss0  = pax[i]*sss0 + Gpx*sss1;
        //const double yss0  = pay[i]*sss0 + Gpy*sss1;
        //const double zss0  = paz[i]*sss0 + Gpz*sss1;
        //const double xss1  = pax[i]*sss1 + Gpx*sss2;
        //const double yss1  = pay[i]*sss1 + Gpy*sss2;
        //const double zss1  = paz[i]*sss1 + Gpz*sss2;
        const double xss2  = pax[i]*sss2 + Gpx*sss3;
        const double yss2  = pay[i]*sss2 + Gpy*sss3;
        const double zss2  = paz[i]*sss2 + Gpz*sss3;
        const double xss3  = pax[i]*sss3 + Gpx*sss4;
        const double yss3  = pay[i]*sss3 + Gpy*sss4;
        const double zss3  = paz[i]*sss3 + Gpz*sss4;
        const double xss4  = pax[i]*sss4 + Gpx*sss5;
        const double yss4  = pay[i]*sss4 + Gpy*sss5;
        const double zss4  = paz[i]*sss4 + Gpz*sss5;
        // [ SS | P ]
        const double ssx1  = Gcx*sss2;
        const double ssy1  = Gcy*sss2;
        const double ssz1  = Gcz*sss2;
        const double ssx2  = Gcx*sss3;
        const double ssy2  = Gcy*sss3;
        const double ssz2  = Gcz*sss3;
        // [ PS | P ]
        const double RGf2  = RzGp5*sss2;
        const double RGf3  = RzGp5*sss3;
        const double xsx1  = Gcx*xss2 + RGf2;
        const double xsy1  = Gcy*xss2;
        const double xsz1  = Gcz*xss2;
        const double ysx1  = Gcx*yss2;
        const double ysy1  = Gcy*yss2 + RGf2;
        const double ysz1  = Gcz*yss2;
        const double zsx1  = Gcx*zss2;
        const double zsy1  = Gcy*zss2;
        const double zsz1  = Gcz*zss2 + RGf2;
        const double xsx2  = Gcx*xss3 + RGf3;
        const double xsy2  = Gcy*xss3;
        const double xsz2  = Gcz*xss3;
        const double ysx2  = Gcx*yss3;
        const double ysy2  = Gcy*yss3 + RGf3;
        const double ysz2  = Gcz*yss3;
        const double zsx2  = Gcx*zss3;
        const double zsy2  = Gcy*zss3;
        const double zsz2  = Gcz*zss3 + RGf3;
        // [ DS | S ]   xy, xz, yz, xx-yy, 3zz-rr
        //const double xxssw0 = pax[i]*xss0 + Gpx*xss1;
        //const double yyssw0 = pay[i]*yss0 + Gpy*yss1;
        //const double xxssw1 = pax[i]*xss1 + Gpx*xss2;
        //const double yyssw1 = pay[i]*yss1 + Gpy*yss2;
        const double xxssw2 = pax[i]*xss2 + Gpx*xss3;
        const double yyssw2 = pay[i]*yss2 + Gpy*yss3;
        const double xxssw3 = pax[i]*xss3 + Gpx*xss4;
        const double yyssw3 = pay[i]*yss3 + Gpy*yss4;
        //const double ass0   = pay[i]*xss0 + Gpy*xss1;
        //const double bss0   = paz[i]*xss0 + Gpz*xss1;
        //const double css0   = paz[i]*yss0 + Gpz*yss1;
        //const double dss0   = 0.5*(xxssw0 - yyssw0);
        //const double ess0   = SQR3I*(paz[i]*zss0 + Gpz*zss1
        //                             - 0.5*(xxssw0 + yyssw0));
        //const double ass1   = pay[i]*xss1 + Gpy*xss2;
        //const double bss1   = paz[i]*xss1 + Gpz*xss2;
        //const double css1   = paz[i]*yss1 + Gpz*yss2;
        //const double dss1   = 0.5*(xxssw1 - yyssw1);
        //const double ess1   = SQR3I*(paz[i]*zss1 + Gpz*zss2
        //                             - 0.5*(xxssw1 + yyssw1));
        const double ass2   = pay[i]*xss2 + Gpy*xss3;
        const double bss2   = paz[i]*xss2 + Gpz*xss3;
        const double css2   = paz[i]*yss2 + Gpz*yss3;
        const double dss2   = 0.5*(xxssw2 - yyssw2);
        const double ess2   = SQR3I*(paz[i]*zss2 + Gpz*zss3
                                     - 0.5*(xxssw2 + yyssw2));
        const double ass3   = pay[i]*xss3 + Gpy*xss4;
        const double bss3   = paz[i]*xss3 + Gpz*xss4;
        const double css3   = paz[i]*yss3 + Gpz*yss4;
        const double dss3   = 0.5*(xxssw3 - yyssw3);
        const double ess3   = SQR3I*(paz[i]*zss3 + Gpz*zss4
                                     - 0.5*(xxssw3 + yyssw3));
        // [ DS | P ]
        const double Rxss2 = RzGp5*xss2;
        const double Ryss2 = RzGp5*yss2;
        const double Rzss2 = RzGp5*zss2;
        const double Rxss3 = RzGp5*xss3;
        const double Ryss3 = RzGp5*yss3;
        const double Rzss3 = RzGp5*zss3;
        const double asx1  = Gcx*ass2 + Ryss2;
        const double asy1  = Gcy*ass2 + Rxss2;
        const double asz1  = Gcz*ass2;
        const double bsx1  = Gcx*bss2 + Rzss2;
        const double bsy1  = Gcy*bss2;
        const double bsz1  = Gcz*bss2 + Rxss2;
        const double csx1  = Gcx*css2;
        const double csy1  = Gcy*css2 + Rzss2;
        const double csz1  = Gcz*css2 + Ryss2;
        const double dsx1  = Gcx*dss2 + Rxss2;
        const double dsy1  = Gcy*dss2 - Ryss2;
        const double dsz1  = Gcz*dss2;
        const double esx1  = Gcx*ess2 - SQR3I *Rxss2;
        const double esy1  = Gcy*ess2 - SQR3I *Ryss2;
        const double esz1  = Gcz*ess2 + SQR3I2*Rzss2;
        const double asx2  = Gcx*ass3 + Ryss3;
        const double asy2  = Gcy*ass3 + Rxss3;
        const double asz2  = Gcz*ass3;
        const double bsx2  = Gcx*bss3 + Rzss3;
        const double bsy2  = Gcy*bss3;
        const double bsz2  = Gcz*bss3 + Rxss3;
        const double csx2  = Gcx*css3;
        const double csy2  = Gcy*css3 + Rzss3;
        const double csz2  = Gcz*css3 + Ryss3;
        const double dsx2  = Gcx*dss3 + Rxss3;
        const double dsy2  = Gcy*dss3 - Ryss3;
        const double dsz2  = Gcz*dss3;
        const double esx2  = Gcx*ess3 - SQR3I *Rxss3;
        const double esy2  = Gcy*ess3 - SQR3I *Ryss3;
        const double esz2  = Gcz*ess3 + SQR3I2*Rzss3;
        // [ PS | D ]
        const double xsx2w1 = Gcx*xsx1 + RzGp5*ssx1;
        const double ysx2w1 = Gcx*ysx1;
        const double zsx2w1 = Gcx*zsx1;
        const double xsy2w1 = Gcy*xsy1;
        const double ysy2w1 = Gcy*ysy1 + RzGp5*ssy1;
        const double zsy2w1 = Gcy*zsy1;
        const double xsz2w1 = Gcz*xsz1;
        const double ysz2w1 = Gcz*ysz1;
        const double zsz2w1 = Gcz*zsz1 + RzGp5*ssz1;
        const double xsx2w2 = Gcx*xsx2 + RzGp5*ssx2;
        const double ysx2w2 = Gcx*ysx2;
        const double zsx2w2 = Gcx*zsx2;
        const double xsy2w2 = Gcy*xsy2;
        const double ysy2w2 = Gcy*ysy2 + RzGp5*ssy2;
        const double zsy2w2 = Gcy*zsy2;
        const double xsz2w2 = Gcz*xsz2;
        const double ysz2w2 = Gcz*ysz2;
        const double zsz2w2 = Gcz*zsz2 + RzGp5*ssz2;

        const double xsa0   = Gcy*xsx1;
        const double ysa0   = Gcx*ysy1;
        const double zsa0   = Gcx*zsy1;
        const double xsb0   = Gcz*xsx1;
        const double ysb0   = Gcx*ysz1;
        const double zsb0   = Gcx*zsz1;
        const double xsc0   = Gcy*xsz1;
        const double ysc0   = Gcz*ysy1;
        const double zsc0   = Gcy*zsz1;
        const double xsd0   = 0.5*(xsx2w1 - xsy2w1);
        const double ysd0   = 0.5*(ysx2w1 - ysy2w1);
        const double zsd0   = 0.5*(zsx2w1 - zsy2w1);
        const double xse0   = SQR3I*(xsz2w1 - 0.5*(xsx2w1 + xsy2w1));
        const double yse0   = SQR3I*(ysz2w1 - 0.5*(ysx2w1 + ysy2w1));
        const double zse0   = SQR3I*(zsz2w1 - 0.5*(zsx2w1 + zsy2w1));
        const double xsa1   = Gcy*xsx2;
        const double ysa1   = Gcx*ysy2;
        const double zsa1   = Gcx*zsy2;
        const double xsb1   = Gcz*xsx2;
        const double ysb1   = Gcx*ysz2;
        const double zsb1   = Gcx*zsz2;
        const double xsc1   = Gcy*xsz2;
        const double ysc1   = Gcz*ysy2;
        const double zsc1   = Gcy*zsz2;
        const double xsd1   = 0.5*(xsx2w2 - xsy2w2);
        const double ysd1   = 0.5*(ysx2w2 - ysy2w2);
        const double zsd1   = 0.5*(zsx2w2 - zsy2w2);
        const double xse1   = SQR3I*(xsz2w2 - 0.5*(xsx2w2 + xsy2w2));
        const double yse1   = SQR3I*(ysz2w2 - 0.5*(ysx2w2 + ysy2w2));
        const double zse1   = SQR3I*(zsz2w2 - 0.5*(zsx2w2 + zsy2w2));
        // [ DS | D ]
        const double Rxsx1  = RzGp5*xsx1;
        const double Rxsy1  = RzGp5*xsy1;
        const double Rxsz1  = RzGp5*xsz1;
        const double Rysx1  = RzGp5*ysx1;
        const double Rysy1  = RzGp5*ysy1;
        const double Rysz1  = RzGp5*ysz1;
        const double Rzsx1  = RzGp5*zsx1;
        const double Rzsy1  = RzGp5*zsy1;
        const double Rzsz1  = RzGp5*zsz1;
        const double Rxsx2  = RzGp5*xsx2;
        const double Rxsy2  = RzGp5*xsy2;
        const double Rxsz2  = RzGp5*xsz2;
        const double Rysx2  = RzGp5*ysx2;
        const double Rysy2  = RzGp5*ysy2;
        const double Rysz2  = RzGp5*ysz2;
        const double Rzsx2  = RzGp5*zsx2;
        const double Rzsy2  = RzGp5*zsy2;
        const double Rzsz2  = RzGp5*zsz2;
        const double asxxw1 = Gcx*asx1 + Rysx1;
        const double asyyw1 = Gcy*asy1 + Rxsy1;
        const double bsxxw1 = Gcx*bsx1 + Rzsx1;
        const double bsyyw1 = Gcy*bsy1;
        const double csxxw1 = Gcx*csx1;
        const double csyyw1 = Gcy*csy1 + Rzsy1;
        const double dsxxw1 = Gcx*dsx1 + Rxsx1;
        const double dsyyw1 = Gcy*dsy1 - Rysy1;
        const double esxxw1 = Gcx*esx1 - SQR3I*Rxsx1;
        const double esyyw1 = Gcy*esy1 - SQR3I*Rysy1;
        const double asxxw2 = Gcx*asx2 + Rysx2;
        const double asyyw2 = Gcy*asy2 + Rxsy2;
        const double bsxxw2 = Gcx*bsx2 + Rzsx2;
        const double bsyyw2 = Gcy*bsy2;
        const double csxxw2 = Gcx*csx2;
        const double csyyw2 = Gcy*csy2 + Rzsy2;
        const double dsxxw2 = Gcx*dsx2 + Rxsx2;
        const double dsyyw2 = Gcy*dsy2 - Rysy2;
        const double esxxw2 = Gcx*esx2 - SQR3I*Rxsx2;
        const double esyyw2 = Gcy*esy2 - SQR3I*Rysy2;

        const double asa0 = Gcx*asy1 + Rysy1;
        const double asb0 = Gcx*asz1 + Rysz1;
        const double asc0 = Gcy*asz1 + Rxsz1;
        const double asd0 = 0.5*(asxxw1 - asyyw1);
        const double ase0 = SQR3I*(Gcz*asz1
                                   - 0.5*(asxxw1 + asyyw1));
        const double bsa0 = Gcx*bsy1 + Rzsy1;
        const double bsb0 = Gcx*bsz1 + Rzsz1;
        const double bsc0 = Gcy*bsz1;
        const double bsd0 = 0.5*(bsxxw1 - bsyyw1);
        const double bse0 = SQR3I*(Gcz*bsz1 + Rxsz1
                                   - 0.5*(bsxxw1 + bsyyw1));
        const double csa0 = Gcx*csy1;
        const double csb0 = Gcx*csz1;
        const double csc0 = Gcy*csz1 + Rzsz1;
        const double csd0 = 0.5*(csxxw1 - csyyw1);
        const double cse0 = SQR3I*(Gcz*csz1 + Rysz1
                                   - 0.5*(csxxw1 + csyyw1));
        const double dsa0 = Gcx*dsy1 + Rxsy1;
        const double dsb0 = Gcx*dsz1 + Rxsz1;
        const double dsc0 = Gcy*dsz1 - Rysz1;
        const double dsd0 = 0.5*(dsxxw1 - dsyyw1);
        const double dse0 = SQR3I*(Gcz*dsz1
                                   - 0.5*(dsxxw1 + dsyyw1));
        const double esa0 = Gcx*esy1 - SQR3I*Rxsy1;
        const double esb0 = Gcx*esz1 - SQR3I*Rxsz1;
        const double esc0 = Gcy*esz1 - SQR3I*Rysz1;
        const double esd0 = 0.5*(esxxw1 - esyyw1);
        const double ese0 = SQR3I*(Gcz*esz1 + SQR3I2*Rzsz1
                                   - 0.5*(esxxw1 + esyyw1));
        const double asa1 = Gcx*asy2 + Rysy2;
        const double asb1 = Gcx*asz2 + Rysz2;
        const double asc1 = Gcy*asz2 + Rxsz2;
        const double asd1 = 0.5*(asxxw2 - asyyw2);
        const double ase1 = SQR3I*(Gcz*asz2
                                   - 0.5*(asxxw2 + asyyw2));
        const double bsa1 = Gcx*bsy2 + Rzsy2;
        const double bsb1 = Gcx*bsz2 + Rzsz2;
        const double bsc1 = Gcy*bsz2;
        const double bsd1 = 0.5*(bsxxw2 - bsyyw2);
        const double bse1 = SQR3I*(Gcz*bsz2 + Rxsz2
                                   - 0.5*(bsxxw2 + bsyyw2));
        const double csa1 = Gcx*csy2;
        const double csb1 = Gcx*csz2;
        const double csc1 = Gcy*csz2 + Rzsz2;
        const double csd1 = 0.5*(csxxw2 - csyyw2);
        const double cse1 = SQR3I*(Gcz*csz2 + Rysz2
                                   - 0.5*(csxxw2 + csyyw2));
        const double dsa1 = Gcx*dsy2 + Rxsy2;
        const double dsb1 = Gcx*dsz2 + Rxsz2;
        const double dsc1 = Gcy*dsz2 - Rysz2;
        const double dsd1 = 0.5*(dsxxw2 - dsyyw2);
        const double dse1 = SQR3I*(Gcz*dsz2
                                   - 0.5*(dsxxw2 + dsyyw2));
        const double esa1 = Gcx*esy2 - SQR3I*Rxsy2;
        const double esb1 = Gcx*esz2 - SQR3I*Rxsz2;
        const double esc1 = Gcy*esz2 - SQR3I*Rysz2;
        const double esd1 = 0.5*(esxxw2 - esyyw2);
        const double ese1 = SQR3I*(Gcz*esz2 + SQR3I2*Rzsz2
                                   - 0.5*(esxxw2 + esyyw2));
        // [ DP | D ]
        const double xsa01 = Zp5*(xsa0 - Rooz*xsa1);
        const double ysa01 = Zp5*(ysa0 - Rooz*ysa1);
        const double zsa01 = Zp5*(zsa0 - Rooz*zsa1);
        const double xsb01 = Zp5*(xsb0 - Rooz*xsb1);
        const double ysb01 = Zp5*(ysb0 - Rooz*ysb1);
        const double zsb01 = Zp5*(zsb0 - Rooz*zsb1);
        const double xsc01 = Zp5*(xsc0 - Rooz*xsc1);
        const double ysc01 = Zp5*(ysc0 - Rooz*ysc1);
        const double zsc01 = Zp5*(zsc0 - Rooz*zsc1);
        const double xsd01 = Zp5*(xsd0 - Rooz*xsd1);
        const double ysd01 = Zp5*(ysd0 - Rooz*ysd1);
        const double zsd01 = Zp5*(zsd0 - Rooz*zsd1);
        const double xse01 = Zp5*(xse0 - Rooz*xse1);
        const double yse01 = Zp5*(yse0 - Rooz*yse1);
        const double zse01 = Zp5*(zse0 - Rooz*zse1);
        const double Rasx1 = RzGp5*asx1;
        const double Rasy1 = RzGp5*asy1;
        const double Rasz1 = RzGp5*asz1;
        const double Rbsx1 = RzGp5*bsx1;
        const double Rbsy1 = RzGp5*bsy1;
        const double Rbsz1 = RzGp5*bsz1;
        const double Rcsx1 = RzGp5*csx1;
        const double Rcsy1 = RzGp5*csy1;
        const double Rcsz1 = RzGp5*csz1;
        const double Rdsx1 = RzGp5*dsx1;
        const double Rdsy1 = RzGp5*dsy1;
        const double Rdsz1 = RzGp5*dsz1;
        const double Resx1 = RzGp5*esx1;
        const double Resy1 = RzGp5*esy1;
        const double Resz1 = RzGp5*esz1;

        ERP[i][ 0] = pbx[i]*asa0 + Gpx*asa1 + ysa01 + Rasy1;
        ERP[i][ 1] = pbx[i]*asb0 + Gpx*asb1 + ysb01 + Rasz1;
        ERP[i][ 2] = pbx[i]*asc0 + Gpx*asc1 + ysc01;
        ERP[i][ 3] = pbx[i]*asd0 + Gpx*asd1 + ysd01 + Rasx1;
        ERP[i][ 4] = pbx[i]*ase0 + Gpx*ase1 + yse01 - SQR3I *Rasx1;
        ERP[i][ 5] = pby[i]*asa0 + Gpy*asa1 + xsa01 + Rasx1;
        ERP[i][ 6] = pby[i]*asb0 + Gpy*asb1 + xsb01;
        ERP[i][ 7] = pby[i]*asc0 + Gpy*asc1 + xsc01 + Rasz1;
        ERP[i][ 8] = pby[i]*asd0 + Gpy*asd1 + xsd01 - Rasy1;
        ERP[i][ 9] = pby[i]*ase0 + Gpy*ase1 + xse01 - SQR3I *Rasy1;
        ERP[i][10] = pbz[i]*asa0 + Gpz*asa1;
        ERP[i][11] = pbz[i]*asb0 + Gpz*asb1         + Rasx1;
        ERP[i][12] = pbz[i]*asc0 + Gpz*asc1         + Rasy1;
        ERP[i][13] = pbz[i]*asd0 + Gpz*asd1;
        ERP[i][14] = pbz[i]*ase0 + Gpz*ase1         + SQR3I2*Rasz1;

        ERP[i][15] = pbx[i]*bsa0 + Gpx*bsa1 + zsa01 + Rbsy1;
        ERP[i][16] = pbx[i]*bsb0 + Gpx*bsb1 + zsb01 + Rbsz1;
        ERP[i][17] = pbx[i]*bsc0 + Gpx*bsc1 + zsc01;
        ERP[i][18] = pbx[i]*bsd0 + Gpx*bsd1 + zsd01 + Rbsx1;
        ERP[i][19] = pbx[i]*bse0 + Gpx*bse1 + zse01 - SQR3I *Rbsx1;
        ERP[i][20] = pby[i]*bsa0 + Gpy*bsa1         + Rbsx1;
        ERP[i][21] = pby[i]*bsb0 + Gpy*bsb1;
        ERP[i][22] = pby[i]*bsc0 + Gpy*bsc1         + Rbsz1;
        ERP[i][23] = pby[i]*bsd0 + Gpy*bsd1         - Rbsy1;
        ERP[i][24] = pby[i]*bse0 + Gpy*bse1         - SQR3I *Rbsy1;
        ERP[i][25] = pbz[i]*bsa0 + Gpz*bsa1 + xsa01;
        ERP[i][26] = pbz[i]*bsb0 + Gpz*bsb1 + xsb01 + Rbsx1;
        ERP[i][27] = pbz[i]*bsc0 + Gpz*bsc1 + xsc01 + Rbsy1;
        ERP[i][28] = pbz[i]*bsd0 + Gpz*bsd1 + xsd01;
        ERP[i][29] = pbz[i]*bse0 + Gpz*bse1 + xse01 + SQR3I2*Rbsz1;

        ERP[i][30] = pbx[i]*csa0 + Gpx*csa1         + Rcsy1;
        ERP[i][31] = pbx[i]*csb0 + Gpx*csb1         + Rcsz1;
        ERP[i][32] = pbx[i]*csc0 + Gpx*csc1;
        ERP[i][33] = pbx[i]*csd0 + Gpx*csd1         + Rcsx1;
        ERP[i][34] = pbx[i]*cse0 + Gpx*cse1         - SQR3I *Rcsx1;
        ERP[i][35] = pby[i]*csa0 + Gpy*csa1 + zsa01 + Rcsx1;
        ERP[i][36] = pby[i]*csb0 + Gpy*csb1 + zsb01;
        ERP[i][37] = pby[i]*csc0 + Gpy*csc1 + zsc01 + Rcsz1;
        ERP[i][38] = pby[i]*csd0 + Gpy*csd1 + zsd01 - Rcsy1;
        ERP[i][39] = pby[i]*cse0 + Gpy*cse1 + zse01 - SQR3I *Rcsy1;
        ERP[i][40] = pbz[i]*csa0 + Gpz*csa1 + ysa01;
        ERP[i][41] = pbz[i]*csb0 + Gpz*csb1 + ysb01 + Rcsx1;
        ERP[i][42] = pbz[i]*csc0 + Gpz*csc1 + ysc01 + Rcsy1;
        ERP[i][43] = pbz[i]*csd0 + Gpz*csd1 + ysd01;
        ERP[i][44] = pbz[i]*cse0 + Gpz*cse1 + yse01 + SQR3I2*Rcsz1;

        ERP[i][45] = pbx[i]*dsa0 + Gpx*dsa1 + xsa01 + Rdsy1;
        ERP[i][46] = pbx[i]*dsb0 + Gpx*dsb1 + xsb01 + Rdsz1;
        ERP[i][47] = pbx[i]*dsc0 + Gpx*dsc1 + xsc01;
        ERP[i][48] = pbx[i]*dsd0 + Gpx*dsd1 + xsd01 + Rdsx1;
        ERP[i][49] = pbx[i]*dse0 + Gpx*dse1 + xse01 - SQR3I *Rdsx1;
        ERP[i][50] = pby[i]*dsa0 + Gpy*dsa1 - ysa01 + Rdsx1;
        ERP[i][51] = pby[i]*dsb0 + Gpy*dsb1 - ysb01;
        ERP[i][52] = pby[i]*dsc0 + Gpy*dsc1 - ysc01 + Rdsz1;
        ERP[i][53] = pby[i]*dsd0 + Gpy*dsd1 - ysd01 - Rdsy1;
        ERP[i][54] = pby[i]*dse0 + Gpy*dse1 - yse01 - SQR3I *Rdsy1;
        ERP[i][55] = pbz[i]*dsa0 + Gpz*dsa1;
        ERP[i][56] = pbz[i]*dsb0 + Gpz*dsb1         + Rdsx1;
        ERP[i][57] = pbz[i]*dsc0 + Gpz*dsc1         + Rdsy1;
        ERP[i][58] = pbz[i]*dsd0 + Gpz*dsd1;
        ERP[i][59] = pbz[i]*dse0 + Gpz*dse1         + SQR3I2*Rdsz1;

        ERP[i][60] = pbx[i]*esa0 + Gpx*esa1 - SQR3I *xsa01 + Resy1;
        ERP[i][61] = pbx[i]*esb0 + Gpx*esb1 - SQR3I *xsb01 + Resz1;
        ERP[i][62] = pbx[i]*esc0 + Gpx*esc1 - SQR3I *xsc01;
        ERP[i][63] = pbx[i]*esd0 + Gpx*esd1 - SQR3I *xsd01 + Resx1;
        ERP[i][64] = pbx[i]*ese0 + Gpx*ese1 - SQR3I *(xse01 + Resx1);
        ERP[i][65] = pby[i]*esa0 + Gpy*esa1 - SQR3I *ysa01 + Resx1;
        ERP[i][66] = pby[i]*esb0 + Gpy*esb1 - SQR3I *ysb01;
        ERP[i][67] = pby[i]*esc0 + Gpy*esc1 - SQR3I *ysc01 + Resz1;
        ERP[i][68] = pby[i]*esd0 + Gpy*esd1 - SQR3I *ysd01 - Resy1;
        ERP[i][69] = pby[i]*ese0 + Gpy*ese1 - SQR3I *(yse01 + Resy1);
        ERP[i][70] = pbz[i]*esa0 + Gpz*esa1 + SQR3I2*zsa01;
        ERP[i][71] = pbz[i]*esb0 + Gpz*esb1 + SQR3I2*zsb01 + Resx1;
        ERP[i][72] = pbz[i]*esc0 + Gpz*esc1 + SQR3I2*zsc01 + Resy1;
        ERP[i][73] = pbz[i]*esd0 + Gpz*esd1 + SQR3I2*zsd01;
        ERP[i][74] = pbz[i]*ese0 + Gpz*ese1 + SQR3I2*(zse01 + Resz1);
    }
}


/*****************************************************************************/
/*                                       */
/*  Primitive 3-center ERI evaluation [ AB | C ] for type [ DD | D ]         */
/*                                       */
/*  np      ; number of first primitive pair data            */
/*  Zpi     ; 1/(Alpha+Beta) from orbital exponents Alpha, Beta  */
/*  px, py, pz  ; (Alpha*ax * Beta*bx) / (Alpha+Beta) , ...          */
/*  pax, pay, paz   ; px-ax, py-ay, pz-az                    */
/*  pbx, pby, pbz   ; px-bx, py-by, pz-bz                    */
/*  HP      ; Zeta**(-3/2) * exp(-(Alpha*Beta/Zeta) * (A-B)**2   */
/*  Gammainv    ; inverse of orbital exponent Gamma for last GTO     */
/*  cx, cy, cz  ; coordinates of last GTO C              */
/*  cc      ; coefficient of last GTO C              */
/*  ERP     ; output primitive ERI                   */
/*          ERP[0,1,2,...,np-1][0] for [ DxyDxy | Dxy ]      */
/*          ERP[0,1,2,...,np-1][1] for [ DxyDxy | Dxz ]      */
/*          ERP[0,1,2,...,np-1][2] for [ DxyDxy | Dyz ]      */
/*          ERP[0,1,2,...,np-1][3] for [ DxyDxy | Dxx-yy ]       */
/*          ERP[0,1,2,...,np-1][4] for [ DxyDxy | D3zz-rr ]      */
/*          ERP[0,1,2,...,np-1][5] for [ DxyDxz | Dxy ]      */
/*          ERP[0,1,2,...,np-1][6] for [ DxyDxz | Dxz ]      */
/*          ERP[0,1,2,...,np-1][7] for [ DxyDxz | Dyz ]      */
/*          ERP[0,1,2,...,np-1][8] for [ DxyDxz | Dxx-yy ]       */
/*          ERP[0,1,2,...,np-1][9] for [ DxyDxz | D3zz-rr ]      */
/*                  .............                */
/*                                       */
/*  FPAI   = 2 * PAI**(5/2)                          */
/*  SQR3I  = 1/sqrt(3)                           */
/*  SQR3I2 = 2/sqrt(3)                           */
/*                                       */
/*****************************************************************************/
void DfEri::eriDDD(const int np,
                   const double* px, const double* py, const double* pz,
                   const double* pax, const double* pay, const double* paz,
                   const double* pbx, const double* pby, const double* pbz,
                   const double cx, const double cy, const double cz,
                   const double cc, const double Gammainv, const double* Zpi, const double* HP, double** ERP)
{
    const double tfw  = TF[6];
    const double rmi0 = RMI[0];
    const double rmi1 = RMI[1];
    const double rmi2 = RMI[2];
    const double rmi3 = RMI[3];
    const double rmi4 = RMI[4];
    const double hc   = FPAI*Gammainv*sqrt(Gammainv)*cc;
    const double Gp5  = 0.5*Gammainv;

    for (int i = 0; i < np; ++i) {
        const double Roo = 1.0 / (Zpi[i] + Gammainv);
        const double pcx = px[i] - cx;
        const double pcy = py[i] - cy;
        const double pcz = pz[i] - cz;
        const double t   = Roo * (pcx*pcx + pcy*pcy + pcz*pcz);

        double f0, f1, f2, f3, f4, f5, f6;
        if (t <= tfw) {
            const int it  = (int)((t+0.015) * D33);
            const double dt  = 0.03*it - t;
            f6  = ((((ADAT[41][it]  *dt + ADAT[40][it]) * dt
                     + ADAT[39][it])*dt + ADAT[38][it]) * dt
                   + ADAT[37][it])*dt + ADAT[36][it];
            const double eed = ((((GA[5]         *dt + GA[4]) * dt
                                  + GA[3])*dt + GA[2]) * dt
                                +   GA[1])*dt + GA[0];
            const double ee  = EDAT[it] * eed;
            const double t2  = 2.0*t;
            f5  = (t2*f6 + ee)*rmi4;
            f4  = (t2*f5 + ee)*rmi3;
            f3  = (t2*f4 + ee)*rmi2;
            f2  = (t2*f3 + ee)*rmi1;
            f1  = (t2*f2 + ee)*rmi0;
            f0  = t2*f1 + ee;
        } else {
            const double tinv = 1.0/t;
            f0   = 0.5*sqrt(M_PI*tinv);
            f1   = 0.5*f0*tinv;
            f2   = 1.5*f1*tinv;
            f3   = 2.5*f2*tinv;
            f4   = 3.5*f3*tinv;
            f5   = 4.5*f4*tinv;
            f6   = 5.5*f5*tinv;
        }

        const double hpq   = hc*sqrt(Roo)*HP[i];
        //const double sss0  = hpq*f0;
        //const double sss1  = hpq*f1;
        const double sss2  = hpq*f2;
        const double sss3  = hpq*f3;
        const double sss4  = hpq*f4;
        const double sss5  = hpq*f5;
        const double sss6  = hpq*f6;
        // [ PS | S ]
        const double RooG  = Roo*Gammainv;
        const double Rooz  = Zpi[i]*Roo;
        const double Zp5   = Zpi[i]*0.5;
        const double RzGp5 = Rooz*Gp5;
        const double Gpx   = -Rooz*pcx;
        const double Gpy   = -Rooz*pcy;
        const double Gpz   = -Rooz*pcz;
        const double Gcx   =  RooG*pcx;
        const double Gcy   =  RooG*pcy;
        const double Gcz   =  RooG*pcz;
        const double abx   = pbx[i] - pax[i];
        const double aby   = pby[i] - pay[i];
        const double abz   = pbz[i] - paz[i];
        const double xss2  = pax[i]*sss2 + Gpx*sss3;
        const double yss2  = pay[i]*sss2 + Gpy*sss3;
        const double zss2  = paz[i]*sss2 + Gpz*sss3;
        const double xss3  = pax[i]*sss3 + Gpx*sss4;
        const double yss3  = pay[i]*sss3 + Gpy*sss4;
        const double zss3  = paz[i]*sss3 + Gpz*sss4;
        const double xss4  = pax[i]*sss4 + Gpx*sss5;
        const double yss4  = pay[i]*sss4 + Gpy*sss5;
        const double zss4  = paz[i]*sss4 + Gpz*sss5;
        const double xss5  = pax[i]*sss5 + Gpx*sss6;
        const double yss5  = pay[i]*sss5 + Gpy*sss6;
        const double zss5  = paz[i]*sss5 + Gpz*sss6;
        // [ SP | S ]
        const double sxs2  = xss2 + abx*sss2;
        const double sys2  = yss2 + aby*sss2;
        const double szs2  = zss2 + abz*sss2;
        // [ PP | S ]
        const double Zps23 = Zp5*(sss2 - Rooz*sss3);
        const double Zps34 = Zp5*(sss3 - Rooz*sss4);
        const double xxs2  = pbx[i]*xss2 + Gpx*xss3 + Zps23;
        const double xys2  = pby[i]*xss2 + Gpy*xss3;
        const double xzs2  = pbz[i]*xss2 + Gpz*xss3;
        const double yxs2  = pbx[i]*yss2 + Gpx*yss3;
        const double yys2  = pby[i]*yss2 + Gpy*yss3 + Zps23;
        const double yzs2  = pbz[i]*yss2 + Gpz*yss3;
        const double zxs2  = pbx[i]*zss2 + Gpx*zss3;
        const double zys2  = pby[i]*zss2 + Gpy*zss3;
        const double zzs2  = pbz[i]*zss2 + Gpz*zss3 + Zps23;
        const double xxs3  = pbx[i]*xss3 + Gpx*xss4 + Zps34;
        const double xys3  = pby[i]*xss3 + Gpy*xss4;
        const double xzs3  = pbz[i]*xss3 + Gpz*xss4;
        const double yxs3  = pbx[i]*yss3 + Gpx*yss4;
        const double yys3  = pby[i]*yss3 + Gpy*yss4 + Zps34;
        const double yzs3  = pbz[i]*yss3 + Gpz*yss4;
        const double zxs3  = pbx[i]*zss3 + Gpx*zss4;
        const double zys3  = pby[i]*zss3 + Gpy*zss4;
        const double zzs3  = pbz[i]*zss3 + Gpz*zss4 + Zps34;
        // [ DS | S ]   xy, xz, yz, xx-yy, 3zz-rr
        const double xxssw2 = pax[i]*xss2 + Gpx*xss3;
        const double yyssw2 = pay[i]*yss2 + Gpy*yss3;
        const double xxssw3 = pax[i]*xss3 + Gpx*xss4;
        const double yyssw3 = pay[i]*yss3 + Gpy*yss4;
        const double xxssw4 = pax[i]*xss4 + Gpx*xss5;
        const double yyssw4 = pay[i]*yss4 + Gpy*yss5;
        const double ass2   = pay[i]*xss2 + Gpy*xss3;
        const double bss2   = paz[i]*xss2 + Gpz*xss3;
        const double css2   = paz[i]*yss2 + Gpz*yss3;
        const double dss2   = 0.5*(xxssw2 - yyssw2);
        const double ess2   = SQR3I*(paz[i]*zss2 + Gpz*zss3
                                     - 0.5*(xxssw2 + yyssw2));
        const double ass3   = pay[i]*xss3 + Gpy*xss4;
        const double bss3   = paz[i]*xss3 + Gpz*xss4;
        const double css3   = paz[i]*yss3 + Gpz*yss4;
        const double dss3   = 0.5*(xxssw3 - yyssw3);
        const double ess3   = SQR3I*(paz[i]*zss3 + Gpz*zss4
                                     - 0.5*(xxssw3 + yyssw3));
        const double ass4   = pay[i]*xss4 + Gpy*xss5;
        const double bss4   = paz[i]*xss4 + Gpz*xss5;
        const double css4   = paz[i]*yss4 + Gpz*yss5;
        const double dss4   = 0.5*(xxssw4 - yyssw4);
        const double ess4   = SQR3I*(paz[i]*zss4 + Gpz*zss5
                                     - 0.5*(xxssw4 + yyssw4));
        // [ DS | S ]   [ DxxS | S ], [ DyyS | S ], [ DzzS | S ]
        const double x2ss2  = pax[i]*xss2 + Gpx*xss3 + Zps23;
        const double y2ss2  = pay[i]*yss2 + Gpy*yss3 + Zps23;
        const double z2ss2  = paz[i]*zss2 + Gpz*zss3 + Zps23;
        const double x2ss3  = pax[i]*xss3 + Gpx*xss4 + Zps34;
        const double y2ss3  = pay[i]*yss3 + Gpy*yss4 + Zps34;
        const double z2ss3  = paz[i]*zss3 + Gpz*zss4 + Zps34;
        // [ SD | S ]
        const double sas2  = xys2 + abx*sys2;
        const double sbs2  = xzs2 + abx*szs2;
        const double scs2  = yzs2 + aby*szs2;
        const double sx2s2 = xxs2 + abx*sxs2;
        const double sy2s2 = yys2 + aby*sys2;
        const double sz2s2 = zzs2 + abz*szs2;
        const double sds2  = 0.5*(sx2s2 - sy2s2);
        const double ses2  = SQR3I*(sz2s2 - 0.5*(sx2s2 + sy2s2));
        // [ DP | S ]
        const double Zxs23 = Zp5*(xss2 - Rooz*xss3);
        const double Zys23 = Zp5*(yss2 - Rooz*yss3);
        const double Zzs23 = Zp5*(zss2 - Rooz*zss3);
        const double Zxs34 = Zp5*(xss3 - Rooz*xss4);
        const double Zys34 = Zp5*(yss3 - Rooz*yss4);
        const double Zzs34 = Zp5*(zss3 - Rooz*zss4);
        const double axs2  = pbx[i]*ass2 + Gpx*ass3 + Zys23;
        const double ays2  = pby[i]*ass2 + Gpy*ass3 + Zxs23;
        const double azs2  = pbz[i]*ass2 + Gpz*ass3;
        const double bxs2  = pbx[i]*bss2 + Gpx*bss3 + Zzs23;
        const double bys2  = pby[i]*bss2 + Gpy*bss3;
        const double bzs2  = pbz[i]*bss2 + Gpz*bss3 + Zxs23;
        const double cxs2  = pbx[i]*css2 + Gpx*css3;
        const double cys2  = pby[i]*css2 + Gpy*css3 + Zzs23;
        const double czs2  = pbz[i]*css2 + Gpz*css3 + Zys23;
        const double dxs2  = pbx[i]*dss2 + Gpx*dss3 + Zxs23;
        const double dys2  = pby[i]*dss2 + Gpy*dss3 - Zys23;
        const double dzs2  = pbz[i]*dss2 + Gpz*dss3;
        const double exs2  = pbx[i]*ess2 + Gpx*ess3 - SQR3I *Zxs23;
        const double eys2  = pby[i]*ess2 + Gpy*ess3 - SQR3I *Zys23;
        const double ezs2  = pbz[i]*ess2 + Gpz*ess3 + SQR3I2*Zzs23;
        const double axs3  = pbx[i]*ass3 + Gpx*ass4 + Zys34;
        const double ays3  = pby[i]*ass3 + Gpy*ass4 + Zxs34;
        const double azs3  = pbz[i]*ass3 + Gpz*ass4;
        const double bxs3  = pbx[i]*bss3 + Gpx*bss4 + Zzs34;
        const double bys3  = pby[i]*bss3 + Gpy*bss4;
        const double bzs3  = pbz[i]*bss3 + Gpz*bss4 + Zxs34;
        const double cxs3  = pbx[i]*css3 + Gpx*css4;
        const double cys3  = pby[i]*css3 + Gpy*css4 + Zzs34;
        const double czs3  = pbz[i]*css3 + Gpz*css4 + Zys34;
        const double dxs3  = pbx[i]*dss3 + Gpx*dss4 + Zxs34;
        const double dys3  = pby[i]*dss3 + Gpy*dss4 - Zys34;
        const double dzs3  = pbz[i]*dss3 + Gpz*dss4;
        const double exs3  = pbx[i]*ess3 + Gpx*ess4 - SQR3I *Zxs34;
        const double eys3  = pby[i]*ess3 + Gpy*ess4 - SQR3I *Zys34;
        const double ezs3  = pbz[i]*ess3 + Gpz*ess4 + SQR3I2*Zzs34;
        // [ DP | S ]   [ DxxPx | S ], [ DyyPy | S ], [ DzzPz | S ]
        const double x2xs2  = pbx[i]*x2ss2 + Gpx*x2ss3 + 2.0*Zxs23;
        const double y2ys2  = pby[i]*y2ss2 + Gpy*y2ss3 + 2.0*Zys23;
        const double z2zs2  = pbz[i]*z2ss2 + Gpz*z2ss3 + 2.0*Zzs23;
        // [ PD | S ]
        const double xas2 = axs2 + aby*xxs2;
        const double xbs2 = bxs2 + abz*xxs2;
        const double xcs2 = azs2 + aby*xzs2;
        const double xds2 = 0.5*(x2xs2 + abx*xxs2 - ays2 - aby*xys2);
        const double xes2 = SQR3I*(bzs2 + abz*xzs2
                                   -0.5*(x2xs2 + abx*xxs2 + ays2 + aby*xys2));
        const double yas2 = ays2 + abx*yys2;
        const double ybs2 = cxs2 + abz*yxs2;
        const double ycs2 = cys2 + abz*yys2;
        const double yds2 = 0.5*(axs2 + abx*yxs2 - y2ys2 - aby*yys2);
        const double yes2 = SQR3I*(czs2 + abz*yzs2
                                   -0.5*(axs2 + abx*yxs2 + y2ys2 + aby*yys2));
        const double zas2 = bys2 + abx*zys2;
        const double zbs2 = bzs2 + abx*zzs2;
        const double zcs2 = czs2 + aby*zzs2;
        const double zds2 = 0.5*(bxs2 + abx*zxs2 - cys2 - aby*zys2);
        const double zes2 = SQR3I*(z2zs2 + abz*zzs2
                                   -0.5*(bxs2 + abx*zxs2 + cys2 + aby*zys2));
        // [ DP | P ]
        const double Rass2 = RzGp5*ass2;
        const double Rbss2 = RzGp5*bss2;
        const double Rcss2 = RzGp5*css2;
        const double Rdss2 = RzGp5*dss2;
        const double Ress2 = RzGp5*ess2;
        const double Rxxs2 = RzGp5*xxs2;
        const double Rxys2 = RzGp5*xys2;
        const double Rxzs2 = RzGp5*xzs2;
        const double Ryxs2 = RzGp5*yxs2;
        const double Ryys2 = RzGp5*yys2;
        const double Ryzs2 = RzGp5*yzs2;
        const double Rzxs2 = RzGp5*zxs2;
        const double Rzys2 = RzGp5*zys2;
        const double Rzzs2 = RzGp5*zzs2;
        // a = xy
        const double axx1 = Gcx*axs2 + Ryxs2 + Rass2;
        const double axy1 = Gcy*axs2 + Rxxs2;
        const double axz1 = Gcz*axs2;
        const double ayx1 = Gcx*ays2 + Ryys2;
        const double ayy1 = Gcy*ays2 + Rxys2 + Rass2;
        const double ayz1 = Gcz*ays2;
        const double azx1 = Gcx*azs2 + Ryzs2;
        const double azy1 = Gcy*azs2 + Rxzs2;
        const double azz1 = Gcz*azs2         + Rass2;
        // b = xz
        const double bxx1 = Gcx*bxs2 + Rzxs2 + Rbss2;
        const double bxy1 = Gcy*bxs2;
        const double bxz1 = Gcz*bxs2 + Rxxs2;
        const double byx1 = Gcx*bys2 + Rzys2;
        const double byy1 = Gcy*bys2         + Rbss2;
        const double byz1 = Gcz*bys2 + Rxys2;
        const double bzx1 = Gcx*bzs2 + Rzzs2;
        const double bzy1 = Gcy*bzs2;
        const double bzz1 = Gcz*bzs2 + Rxzs2 + Rbss2;
        // c = yz
        const double cxx1 = Gcx*cxs2         + Rcss2;
        const double cxy1 = Gcy*cxs2 + Rzxs2;
        const double cxz1 = Gcz*cxs2 + Ryxs2;
        const double cyx1 = Gcx*cys2;
        const double cyy1 = Gcy*cys2 + Rzys2 + Rcss2;
        const double cyz1 = Gcz*cys2 + Ryys2;
        const double czx1 = Gcx*czs2;
        const double czy1 = Gcy*czs2 + Rzzs2;
        const double czz1 = Gcz*czs2 + Ryzs2 + Rcss2;
        // d = xx-yy
        const double dxx1 = Gcx*dxs2 + Rxxs2 + Rdss2;
        const double dxy1 = Gcy*dxs2 - Ryxs2;
        const double dxz1 = Gcz*dxs2;
        const double dyx1 = Gcx*dys2 + Rxys2;
        const double dyy1 = Gcy*dys2 - Ryys2 + Rdss2;
        const double dyz1 = Gcz*dys2;
        const double dzx1 = Gcx*dzs2 + Rxzs2;
        const double dzy1 = Gcy*dzs2 - Ryzs2;
        const double dzz1 = Gcz*dzs2         + Rdss2;
        // e = 3zz-rr
        const double exx1 = Gcx*exs2 - SQR3I *Rxxs2 + Ress2;
        const double exy1 = Gcy*exs2 - SQR3I *Ryxs2;
        const double exz1 = Gcz*exs2 + SQR3I2*Rzxs2;
        const double eyx1 = Gcx*eys2 - SQR3I *Rxys2;
        const double eyy1 = Gcy*eys2 - SQR3I *Ryys2 + Ress2;
        const double eyz1 = Gcz*eys2 + SQR3I2*Rzys2;
        const double ezx1 = Gcx*ezs2 - SQR3I *Rxzs2;
        const double ezy1 = Gcy*ezs2 - SQR3I *Ryzs2;
        const double ezz1 = Gcz*ezs2 + SQR3I2*Rzzs2 + Ress2;
        // [ PD | P ]
        const double Rsas2 = RzGp5*sas2;
        const double Rsbs2 = RzGp5*sbs2;
        const double Rscs2 = RzGp5*scs2;
        const double Rsds2 = RzGp5*sds2;
        const double Rses2 = RzGp5*ses2;
        //const double Rxxs2 = RzGp5*xxs2;
        //const double Rxys2 = RzGp5*xys2;
        //const double Rxzs2 = RzGp5*xzs2;
        //const double Ryxs2 = RzGp5*yxs2;
        //const double Ryys2 = RzGp5*yys2;
        //const double Ryzs2 = RzGp5*yzs2;
        //const double Rzxs2 = RzGp5*zxs2;
        //const double Rzys2 = RzGp5*zys2;
        //const double Rzzs2 = RzGp5*zzs2;
        const double xax1  = Gcx*xas2 + Rsas2 + Rxys2;
        const double xay1  = Gcy*xas2         + Rxxs2;
        const double xaz1  = Gcz*xas2;
        const double xbx1  = Gcx*xbs2 + Rsbs2 + Rxzs2;
        const double xby1  = Gcy*xbs2;
        const double xbz1  = Gcz*xbs2         + Rxxs2;
        const double xcx1  = Gcx*xcs2 + Rscs2;
        const double xcy1  = Gcy*xcs2         + Rxzs2;
        const double xcz1  = Gcz*xcs2         + Rxys2;
        const double xdx1  = Gcx*xds2 + Rsds2 + Rxxs2;
        const double xdy1  = Gcy*xds2         - Rxys2;
        const double xdz1  = Gcz*xds2;
        const double xex1  = Gcx*xes2 + Rses2 - SQR3I *Rxxs2;
        const double xey1  = Gcy*xes2         - SQR3I *Rxys2;
        const double xez1  = Gcz*xes2         + SQR3I2*Rxzs2;
        const double yax1  = Gcx*yas2         + Ryys2;
        const double yay1  = Gcy*yas2 + Rsas2 + Ryxs2;
        const double yaz1  = Gcz*yas2;
        const double ybx1  = Gcx*ybs2         + Ryzs2;
        const double yby1  = Gcy*ybs2 + Rsbs2;
        const double ybz1  = Gcz*ybs2         + Ryxs2;
        const double ycx1  = Gcx*ycs2;
        const double ycy1  = Gcy*ycs2 + Rscs2 + Ryzs2;
        const double ycz1  = Gcz*ycs2         + Ryys2;
        const double ydx1  = Gcx*yds2         + Ryxs2;
        const double ydy1  = Gcy*yds2 + Rsds2 - Ryys2;
        const double ydz1  = Gcz*yds2;
        const double yex1  = Gcx*yes2         - SQR3I *Ryxs2;
        const double yey1  = Gcy*yes2 + Rses2 - SQR3I *Ryys2;
        const double yez1  = Gcz*yes2         + SQR3I2*Ryzs2;
        const double zax1  = Gcx*zas2         + Rzys2;
        const double zay1  = Gcy*zas2         + Rzxs2;
        const double zaz1  = Gcz*zas2 + Rsas2;
        const double zbx1  = Gcx*zbs2         + Rzzs2;
        const double zby1  = Gcy*zbs2;
        const double zbz1  = Gcz*zbs2 + Rsbs2 + Rzxs2;
        const double zcx1  = Gcx*zcs2;
        const double zcy1  = Gcy*zcs2         + Rzzs2;
        const double zcz1  = Gcz*zcs2 + Rscs2 + Rzys2;
        const double zdx1  = Gcx*zds2         + Rzxs2;
        const double zdy1  = Gcy*zds2         - Rzys2;
        const double zdz1  = Gcz*zds2 + Rsds2;
        const double zex1  = Gcx*zes2         - SQR3I *Rzxs2;
        const double zey1  = Gcy*zes2         - SQR3I *Rzys2;
        const double zez1  = Gcz*zes2 + Rses2 + SQR3I2*Rzzs2;
        // [ DD | S ]
        const double xxs23 = Zp5*(xxs2 - Rooz*xxs3);
        const double xys23 = Zp5*(xys2 - Rooz*xys3);
        const double xzs23 = Zp5*(xzs2 - Rooz*xzs3);
        const double yxs23 = Zp5*(yxs2 - Rooz*yxs3);
        const double yys23 = Zp5*(yys2 - Rooz*yys3);
        const double yzs23 = Zp5*(yzs2 - Rooz*yzs3);
        const double zxs23 = Zp5*(zxs2 - Rooz*zxs3);
        const double zys23 = Zp5*(zys2 - Rooz*zys3);
        const double zzs23 = Zp5*(zzs2 - Rooz*zzs3);
        const double axxsw = pbx[i]*axs2 + Gpx*axs3 + yxs23;
        const double ayysw = pby[i]*ays2 + Gpy*ays3 + xys23;
        const double azzsw = pbz[i]*azs2 + Gpz*azs3;
        const double bxxsw = pbx[i]*bxs2 + Gpx*bxs3 + zxs23;
        const double byysw = pby[i]*bys2 + Gpy*bys3;
        const double bzzsw = pbz[i]*bzs2 + Gpz*bzs3 + xzs23;
        const double cxxsw = pbx[i]*cxs2 + Gpx*cxs3;
        const double cyysw = pby[i]*cys2 + Gpy*cys3 + zys23;
        const double czzsw = pbz[i]*czs2 + Gpz*czs3 + yzs23;
        const double dxxsw = pbx[i]*dxs2 + Gpx*dxs3 + xxs23;
        const double dyysw = pby[i]*dys2 + Gpy*dys3 - yys23;
        const double dzzsw = pbz[i]*dzs2 + Gpz*dzs3;
        const double exxsw = pbx[i]*exs2 + Gpx*exs3 - SQR3I *xxs23;
        const double eyysw = pby[i]*eys2 + Gpy*eys3 - SQR3I *yys23;
        const double ezzsw = pbz[i]*ezs2 + Gpz*ezs3 + SQR3I2*zzs23;

        const double aas2  = pbx[i]*ays2 + Gpx*ays3 + yys23;
        const double abs2  = pbx[i]*azs2 + Gpx*azs3 + yzs23;
        const double acs2  = pby[i]*azs2 + Gpy*azs3 + xzs23;
        const double ads2  = 0.5*(axxsw - ayysw);
        const double aes2  = SQR3I*(azzsw - 0.5*(axxsw + ayysw));
        const double bas2  = pbx[i]*bys2 + Gpx*bys3 + zys23;
        const double bbs2  = pbx[i]*bzs2 + Gpx*bzs3 + zzs23;
        const double bcs2  = pby[i]*bzs2 + Gpy*bzs3;
        const double bds2  = 0.5*(bxxsw - byysw);
        const double bes2  = SQR3I*(bzzsw - 0.5*(bxxsw + byysw));
        const double cas2  = pbx[i]*cys2 + Gpx*cys3;
        const double cbs2  = pbx[i]*czs2 + Gpx*czs3;
        const double ccs2  = pby[i]*czs2 + Gpy*czs3 + zzs23;
        const double cds2  = 0.5*(cxxsw - cyysw);
        const double ces2  = SQR3I*(czzsw - 0.5*(cxxsw + cyysw));
        const double das2  = pbx[i]*dys2 + Gpx*dys3 + xys23;
        const double dbs2  = pbx[i]*dzs2 + Gpx*dzs3 + xzs23;
        const double dcs2  = pby[i]*dzs2 + Gpy*dzs3 - yzs23;
        const double dds2  = 0.5*(dxxsw - dyysw);
        const double des2  = SQR3I*(dzzsw - 0.5*(dxxsw + dyysw));
        const double eas2  = pbx[i]*eys2 + Gpx*eys3 - SQR3I*xys23;
        const double ebs2  = pbx[i]*ezs2 + Gpx*ezs3 - SQR3I*xzs23;
        const double ecs2  = pby[i]*ezs2 + Gpy*ezs3 - SQR3I*yzs23;
        const double eds2  = 0.5*(exxsw - eyysw);
        const double ees2  = SQR3I*(ezzsw - 0.5*(exxsw + eyysw));
        // [ DD | P ]
        const double Raxs2 = RzGp5*axs2;
        const double Rays2 = RzGp5*ays2;
        const double Razs2 = RzGp5*azs2;
        const double Rbxs2 = RzGp5*bxs2;
        const double Rbys2 = RzGp5*bys2;
        const double Rbzs2 = RzGp5*bzs2;
        const double Rcxs2 = RzGp5*cxs2;
        const double Rcys2 = RzGp5*cys2;
        const double Rczs2 = RzGp5*czs2;
        const double Rdxs2 = RzGp5*dxs2;
        const double Rdys2 = RzGp5*dys2;
        const double Rdzs2 = RzGp5*dzs2;
        const double Rexs2 = RzGp5*exs2;
        const double Reys2 = RzGp5*eys2;
        const double Rezs2 = RzGp5*ezs2;
        const double Rxas2 = RzGp5*xas2;
        const double Rxbs2 = RzGp5*xbs2;
        const double Rxcs2 = RzGp5*xcs2;
        const double Rxds2 = RzGp5*xds2;
        const double Rxes2 = RzGp5*xes2;
        const double Ryas2 = RzGp5*yas2;
        const double Rybs2 = RzGp5*ybs2;
        const double Rycs2 = RzGp5*ycs2;
        const double Ryds2 = RzGp5*yds2;
        const double Ryes2 = RzGp5*yes2;
        const double Rzas2 = RzGp5*zas2;
        const double Rzbs2 = RzGp5*zbs2;
        const double Rzcs2 = RzGp5*zcs2;
        const double Rzds2 = RzGp5*zds2;
        const double Rzes2 = RzGp5*zes2;

        const double aax1 = Gcx*aas2 + Ryas2 + Rays2;
        const double aay1 = Gcy*aas2 + Rxas2 + Raxs2;
        const double aaz1 = Gcz*aas2;
        const double abx1 = Gcx*abs2 + Rybs2 + Razs2;
        const double aby1 = Gcy*abs2 + Rxbs2;
        const double abz1 = Gcz*abs2         + Raxs2;
        const double acx1 = Gcx*acs2 + Rycs2;
        const double acy1 = Gcy*acs2 + Rxcs2 + Razs2;
        const double acz1 = Gcz*acs2         + Rays2;
        const double adx1 = Gcx*ads2 + Ryds2 + Raxs2;
        const double ady1 = Gcy*ads2 + Rxds2 - Rays2;
        const double adz1 = Gcz*ads2;
        const double aex1 = Gcx*aes2 + Ryes2 - SQR3I *Raxs2;
        const double aey1 = Gcy*aes2 + Rxes2 - SQR3I *Rays2;
        const double aez1 = Gcz*aes2         + SQR3I2*Razs2;

        const double bax1 = Gcx*bas2 + Rzas2 + Rbys2;
        const double bay1 = Gcy*bas2         + Rbxs2;
        const double baz1 = Gcz*bas2 + Rxas2;
        const double bbx1 = Gcx*bbs2 + Rzbs2 + Rbzs2;
        const double bby1 = Gcy*bbs2;
        const double bbz1 = Gcz*bbs2 + Rxbs2 + Rbxs2;
        const double bcx1 = Gcx*bcs2 + Rzcs2;
        const double bcy1 = Gcy*bcs2         + Rbzs2;
        const double bcz1 = Gcz*bcs2 + Rxcs2 + Rbys2;
        const double bdx1 = Gcx*bds2 + Rzds2 + Rbxs2;
        const double bdy1 = Gcy*bds2         - Rbys2;
        const double bdz1 = Gcz*bds2 + Rxds2;
        const double bex1 = Gcx*bes2 + Rzes2 - SQR3I *Rbxs2;
        const double bey1 = Gcy*bes2         - SQR3I *Rbys2;
        const double bez1 = Gcz*bes2 + Rxes2 + SQR3I2*Rbzs2;

        const double cax1 = Gcx*cas2         + Rcys2;
        const double cay1 = Gcy*cas2 + Rzas2 + Rcxs2;
        const double caz1 = Gcz*cas2 + Ryas2;
        const double cbx1 = Gcx*cbs2         + Rczs2;
        const double cby1 = Gcy*cbs2 + Rzbs2;
        const double cbz1 = Gcz*cbs2 + Rybs2 + Rcxs2;
        const double ccx1 = Gcx*ccs2;
        const double ccy1 = Gcy*ccs2 + Rzcs2 + Rczs2;
        const double ccz1 = Gcz*ccs2 + Rycs2 + Rcys2;
        const double cdx1 = Gcx*cds2         + Rcxs2;
        const double cdy1 = Gcy*cds2 + Rzds2 - Rcys2;
        const double cdz1 = Gcz*cds2 + Ryds2;
        const double cex1 = Gcx*ces2         - SQR3I *Rcxs2;
        const double cey1 = Gcy*ces2 + Rzes2 - SQR3I *Rcys2;
        const double cez1 = Gcz*ces2 + Ryes2 + SQR3I2*Rczs2;

        const double dax1 = Gcx*das2 + Rxas2 + Rdys2;
        const double day1 = Gcy*das2 - Ryas2 + Rdxs2;
        const double daz1 = Gcz*das2;
        const double dbx1 = Gcx*dbs2 + Rxbs2 + Rdzs2;
        const double dby1 = Gcy*dbs2 - Rybs2;
        const double dbz1 = Gcz*dbs2         + Rdxs2;
        const double dcx1 = Gcx*dcs2 + Rxcs2;
        const double dcy1 = Gcy*dcs2 - Rycs2 + Rdzs2;
        const double dcz1 = Gcz*dcs2         + Rdys2;
        const double ddx1 = Gcx*dds2 + Rxds2 + Rdxs2;
        const double ddy1 = Gcy*dds2 - Ryds2 - Rdys2;
        const double ddz1 = Gcz*dds2;
        const double dex1 = Gcx*des2 + Rxes2 - SQR3I *Rdxs2;
        const double dey1 = Gcy*des2 - Ryes2 - SQR3I *Rdys2;
        const double dez1 = Gcz*des2         + SQR3I2*Rdzs2;

        const double eax1 = Gcx*eas2 - SQR3I *Rxas2 + Reys2;
        const double eay1 = Gcy*eas2 - SQR3I *Ryas2 + Rexs2;
        const double eaz1 = Gcz*eas2 + SQR3I2*Rzas2;
        const double ebx1 = Gcx*ebs2 - SQR3I *Rxbs2 + Rezs2;
        const double eby1 = Gcy*ebs2 - SQR3I *Rybs2;
        const double ebz1 = Gcz*ebs2 + SQR3I2*Rzbs2 + Rexs2;
        const double ecx1 = Gcx*ecs2 - SQR3I *Rxcs2;
        const double ecy1 = Gcy*ecs2 - SQR3I *Rycs2 + Rezs2;
        const double ecz1 = Gcz*ecs2 + SQR3I2*Rzcs2 + Reys2;
        const double edx1 = Gcx*eds2 - SQR3I *Rxds2 + Rexs2;
        const double edy1 = Gcy*eds2 - SQR3I *Ryds2 - Reys2;
        const double edz1 = Gcz*eds2 + SQR3I2*Rzds2;
        const double eex1 = Gcx*ees2 - SQR3I *(Rxes2 + Rexs2);
        const double eey1 = Gcy*ees2 - SQR3I *(Ryes2 + Reys2);
        const double eez1 = Gcz*ees2 + SQR3I2*(Rzes2 + Rezs2);
        // [ DD | D ]
        const double Rxax1 = RzGp5*xax1;
        const double Rxay1 = RzGp5*xay1;
        const double Rxaz1 = RzGp5*xaz1;
        const double Rxbx1 = RzGp5*xbx1;
        const double Rxby1 = RzGp5*xby1;
        const double Rxbz1 = RzGp5*xbz1;
        const double Rxcx1 = RzGp5*xcx1;
        const double Rxcy1 = RzGp5*xcy1;
        const double Rxcz1 = RzGp5*xcz1;
        const double Rxdx1 = RzGp5*xdx1;
        const double Rxdy1 = RzGp5*xdy1;
        const double Rxdz1 = RzGp5*xdz1;
        const double Rxex1 = RzGp5*xex1;
        const double Rxey1 = RzGp5*xey1;
        const double Rxez1 = RzGp5*xez1;
        const double Ryax1 = RzGp5*yax1;
        const double Ryay1 = RzGp5*yay1;
        const double Ryaz1 = RzGp5*yaz1;
        const double Rybx1 = RzGp5*ybx1;
        const double Ryby1 = RzGp5*yby1;
        const double Rybz1 = RzGp5*ybz1;
        const double Rycx1 = RzGp5*ycx1;
        const double Rycy1 = RzGp5*ycy1;
        const double Rycz1 = RzGp5*ycz1;
        const double Rydx1 = RzGp5*ydx1;
        const double Rydy1 = RzGp5*ydy1;
        const double Rydz1 = RzGp5*ydz1;
        const double Ryex1 = RzGp5*yex1;
        const double Ryey1 = RzGp5*yey1;
        const double Ryez1 = RzGp5*yez1;
        const double Rzax1 = RzGp5*zax1;
        const double Rzay1 = RzGp5*zay1;
        const double Rzaz1 = RzGp5*zaz1;
        const double Rzbx1 = RzGp5*zbx1;
        const double Rzby1 = RzGp5*zby1;
        const double Rzbz1 = RzGp5*zbz1;
        const double Rzcx1 = RzGp5*zcx1;
        const double Rzcy1 = RzGp5*zcy1;
        const double Rzcz1 = RzGp5*zcz1;
        const double Rzdx1 = RzGp5*zdx1;
        const double Rzdy1 = RzGp5*zdy1;
        const double Rzdz1 = RzGp5*zdz1;
        const double Rzex1 = RzGp5*zex1;
        const double Rzey1 = RzGp5*zey1;
        const double Rzez1 = RzGp5*zez1;

        const double Raxx1 = RzGp5*axx1;
        const double Raxy1 = RzGp5*axy1;
        const double Raxz1 = RzGp5*axz1;
        const double Rayx1 = RzGp5*ayx1;
        const double Rayy1 = RzGp5*ayy1;
        const double Rayz1 = RzGp5*ayz1;
        const double Razx1 = RzGp5*azx1;
        const double Razy1 = RzGp5*azy1;
        const double Razz1 = RzGp5*azz1;
        const double Rbxx1 = RzGp5*bxx1;
        const double Rbxy1 = RzGp5*bxy1;
        const double Rbxz1 = RzGp5*bxz1;
        const double Rbyx1 = RzGp5*byx1;
        const double Rbyy1 = RzGp5*byy1;
        const double Rbyz1 = RzGp5*byz1;
        const double Rbzx1 = RzGp5*bzx1;
        const double Rbzy1 = RzGp5*bzy1;
        const double Rbzz1 = RzGp5*bzz1;
        const double Rcxx1 = RzGp5*cxx1;
        const double Rcxy1 = RzGp5*cxy1;
        const double Rcxz1 = RzGp5*cxz1;
        const double Rcyx1 = RzGp5*cyx1;
        const double Rcyy1 = RzGp5*cyy1;
        const double Rcyz1 = RzGp5*cyz1;
        const double Rczx1 = RzGp5*czx1;
        const double Rczy1 = RzGp5*czy1;
        const double Rczz1 = RzGp5*czz1;
        const double Rdxx1 = RzGp5*dxx1;
        const double Rdxy1 = RzGp5*dxy1;
        const double Rdxz1 = RzGp5*dxz1;
        const double Rdyx1 = RzGp5*dyx1;
        const double Rdyy1 = RzGp5*dyy1;
        const double Rdyz1 = RzGp5*dyz1;
        const double Rdzx1 = RzGp5*dzx1;
        const double Rdzy1 = RzGp5*dzy1;
        const double Rdzz1 = RzGp5*dzz1;
        const double Rexx1 = RzGp5*exx1;
        const double Rexy1 = RzGp5*exy1;
        const double Rexz1 = RzGp5*exz1;
        const double Reyx1 = RzGp5*eyx1;
        const double Reyy1 = RzGp5*eyy1;
        const double Reyz1 = RzGp5*eyz1;
        const double Rezx1 = RzGp5*ezx1;
        const double Rezy1 = RzGp5*ezy1;
        const double Rezz1 = RzGp5*ezz1;

        const double aaxxw = Gcx*aax1 + Ryax1 + Rayx1;
        const double aayyw = Gcy*aay1 + Rxay1 + Raxy1;
        const double aazzw = Gcz*aaz1;
        const double abxxw = Gcx*abx1 + Rybx1 + Razx1;
        const double abyyw = Gcy*aby1 + Rxby1;
        const double abzzw = Gcz*abz1         + Raxz1;
        const double acxxw = Gcx*acx1 + Rycx1;
        const double acyyw = Gcy*acy1 + Rxcy1 + Razy1;
        const double aczzw = Gcz*acz1         + Rayz1;
        const double adxxw = Gcx*adx1 + Rydx1 + Raxx1;
        const double adyyw = Gcy*ady1 + Rxdy1 - Rayy1;
        const double adzzw = Gcz*adz1;
        const double aexxw = Gcx*aex1 + Ryex1 - SQR3I *Raxx1;
        const double aeyyw = Gcy*aey1 + Rxey1 - SQR3I *Rayy1;
        const double aezzw = Gcz*aez1         + SQR3I2*Razz1;

        const double baxxw = Gcx*bax1 + Rzax1 + Rbyx1;
        const double bayyw = Gcy*bay1         + Rbxy1;
        const double bazzw = Gcz*baz1 + Rxaz1;
        const double bbxxw = Gcx*bbx1 + Rzbx1 + Rbzx1;
        const double bbyyw = Gcy*bby1;
        const double bbzzw = Gcz*bbz1 + Rxbz1 + Rbxz1;
        const double bcxxw = Gcx*bcx1 + Rzcx1;
        const double bcyyw = Gcy*bcy1         + Rbzy1;
        const double bczzw = Gcz*bcz1 + Rxcz1 + Rbyz1;
        const double bdxxw = Gcx*bdx1 + Rzdx1 + Rbxx1;
        const double bdyyw = Gcy*bdy1         - Rbyy1;
        const double bdzzw = Gcz*bdz1 + Rxdz1;
        const double bexxw = Gcx*bex1 + Rzex1 - SQR3I *Rbxx1;
        const double beyyw = Gcy*bey1         - SQR3I *Rbyy1;
        const double bezzw = Gcz*bez1 + Rxez1 + SQR3I2*Rbzz1;

        const double caxxw = Gcx*cax1         + Rcyx1;
        const double cayyw = Gcy*cay1 + Rzay1 + Rcxy1;
        const double cazzw = Gcz*caz1 + Ryaz1;
        const double cbxxw = Gcx*cbx1         + Rczx1;
        const double cbyyw = Gcy*cby1 + Rzby1;
        const double cbzzw = Gcz*cbz1 + Rybz1 + Rcxz1;
        const double ccxxw = Gcx*ccx1;
        const double ccyyw = Gcy*ccy1 + Rzcy1 + Rczy1;
        const double cczzw = Gcz*ccz1 + Rycz1 + Rcyz1;
        const double cdxxw = Gcx*cdx1         + Rcxx1;
        const double cdyyw = Gcy*cdy1 + Rzdy1 - Rcyy1;
        const double cdzzw = Gcz*cdz1 + Rydz1;
        const double cexxw = Gcx*cex1         - SQR3I *Rcxx1;
        const double ceyyw = Gcy*cey1 + Rzey1 - SQR3I *Rcyy1;
        const double cezzw = Gcz*cez1 + Ryez1 + SQR3I2*Rczz1;

        const double daxxw = Gcx*dax1 + Rxax1 + Rdyx1;
        const double dayyw = Gcy*day1 - Ryay1 + Rdxy1;
        const double dazzw = Gcz*daz1;
        const double dbxxw = Gcx*dbx1 + Rxbx1 + Rdzx1;
        const double dbyyw = Gcy*dby1 - Ryby1;
        const double dbzzw = Gcz*dbz1         + Rdxz1;
        const double dcxxw = Gcx*dcx1 + Rxcx1;
        const double dcyyw = Gcy*dcy1 - Rycy1 + Rdzy1;
        const double dczzw = Gcz*dcz1         + Rdyz1;
        const double ddxxw = Gcx*ddx1 + Rxdx1 + Rdxx1;
        const double ddyyw = Gcy*ddy1 - Rydy1 - Rdyy1;
        const double ddzzw = Gcz*ddz1;
        const double dexxw = Gcx*dex1 + Rxex1 - SQR3I *Rdxx1;
        const double deyyw = Gcy*dey1 - Ryey1 - SQR3I *Rdyy1;
        const double dezzw = Gcz*dez1         + SQR3I2*Rdzz1;

        const double eaxxw = Gcx*eax1 - SQR3I *Rxax1 + Reyx1;
        const double eayyw = Gcy*eay1 - SQR3I *Ryay1 + Rexy1;
        const double eazzw = Gcz*eaz1 + SQR3I2*Rzaz1;
        const double ebxxw = Gcx*ebx1 - SQR3I *Rxbx1 + Rezx1;
        const double ebyyw = Gcy*eby1 - SQR3I *Ryby1;
        const double ebzzw = Gcz*ebz1 + SQR3I2*Rzbz1 + Rexz1;
        const double ecxxw = Gcx*ecx1 - SQR3I *Rxcx1;
        const double ecyyw = Gcy*ecy1 - SQR3I *Rycy1 + Rezy1;
        const double eczzw = Gcz*ecz1 + SQR3I2*Rzcz1 + Reyz1;
        const double edxxw = Gcx*edx1 - SQR3I *Rxdx1 + Rexx1;
        const double edyyw = Gcy*edy1 - SQR3I *Rydy1 - Reyy1;
        const double edzzw = Gcz*edz1 + SQR3I2*Rzdz1;
        const double eexxw = Gcx*eex1 - SQR3I *(Rxex1 + Rexx1);
        const double eeyyw = Gcy*eey1 - SQR3I *(Ryey1 + Reyy1);
        const double eezzw = Gcz*eez1 + SQR3I2*(Rzez1 + Rezz1);

        ERP[i][  0] = Gcx*aay1 + Ryay1 + Rayy1;
        ERP[i][  1] = Gcz*aax1;
        ERP[i][  2] = Gcz*aay1;
        ERP[i][  3] = 0.5*(aaxxw - aayyw);
        ERP[i][  4] = SQR3I*(aazzw - 0.5*(aaxxw + aayyw));
        ERP[i][  5] = Gcx*aby1 + Ryby1 + Razy1;
        ERP[i][  6] = Gcz*abx1         + Raxx1;
        ERP[i][  7] = Gcz*aby1         + Raxy1;
        ERP[i][  8] = 0.5*(abxxw - abyyw);
        ERP[i][  9] = SQR3I*(abzzw - 0.5*(abxxw + abyyw));
        ERP[i][ 10] = Gcx*acy1 + Rycy1;
        ERP[i][ 11] = Gcz*acx1         + Rayx1;
        ERP[i][ 12] = Gcz*acy1         + Rayy1;
        ERP[i][ 13] = 0.5*(acxxw - acyyw);
        ERP[i][ 14] = SQR3I*(aczzw - 0.5*(acxxw + acyyw));
        ERP[i][ 15] = Gcx*ady1 + Rydy1 + Raxy1;
        ERP[i][ 16] = Gcz*adx1;
        ERP[i][ 17] = Gcz*ady1;
        ERP[i][ 18] = 0.5*(adxxw - adyyw);
        ERP[i][ 19] = SQR3I*(adzzw - 0.5*(adxxw + adyyw));
        ERP[i][ 20] = Gcx*aey1 + Ryey1 - SQR3I *Raxy1;
        ERP[i][ 21] = Gcz*aex1         + SQR3I2*Razx1;
        ERP[i][ 22] = Gcz*aey1         + SQR3I2*Razy1;
        ERP[i][ 23] = 0.5*(aexxw - aeyyw);
        ERP[i][ 24] = SQR3I*(aezzw - 0.5*(aexxw + aeyyw));

        ERP[i][ 25] = Gcx*bay1 + Rzay1 + Rbyy1;
        ERP[i][ 26] = Gcz*bax1 + Rxax1;
        ERP[i][ 27] = Gcz*bay1 + Rxay1;
        ERP[i][ 28] = 0.5*(baxxw - bayyw);
        ERP[i][ 29] = SQR3I*(bazzw - 0.5*(baxxw + bayyw));
        ERP[i][ 30] = Gcx*bby1 + Rzby1 + Rbzy1;
        ERP[i][ 31] = Gcz*bbx1 + Rxbx1 + Rbxx1;
        ERP[i][ 32] = Gcz*bby1 + Rxby1 + Rbxy1;
        ERP[i][ 33] = 0.5*(bbxxw - bbyyw);
        ERP[i][ 34] = SQR3I*(bbzzw - 0.5*(bbxxw + bbyyw));
        ERP[i][ 35] = Gcx*bcy1 + Rzcy1;
        ERP[i][ 36] = Gcz*bcx1 + Rxcx1 + Rbyx1;
        ERP[i][ 37] = Gcz*bcy1 + Rxcy1 + Rbyy1;
        ERP[i][ 38] = 0.5*(bcxxw - bcyyw);
        ERP[i][ 39] = SQR3I*(bczzw - 0.5*(bcxxw + bcyyw));
        ERP[i][ 40] = Gcx*bdy1 + Rzdy1 + Rbxy1;
        ERP[i][ 41] = Gcz*bdx1 + Rxdx1;
        ERP[i][ 42] = Gcz*bdy1 + Rxdy1;
        ERP[i][ 43] = 0.5*(bdxxw - bdyyw);
        ERP[i][ 44] = SQR3I*(bdzzw - 0.5*(bdxxw + bdyyw));
        ERP[i][ 45] = Gcx*bey1 + Rzey1 - SQR3I *Rbxy1;
        ERP[i][ 46] = Gcz*bex1 + Rxex1 + SQR3I2*Rbzx1;
        ERP[i][ 47] = Gcz*bey1 + Rxey1 + SQR3I2*Rbzy1;
        ERP[i][ 48] = 0.5*(bexxw - beyyw);
        ERP[i][ 49] = SQR3I*(bezzw - 0.5*(bexxw + beyyw));

        ERP[i][ 50] = Gcx*cay1         + Rcyy1;
        ERP[i][ 51] = Gcz*cax1 + Ryax1;
        ERP[i][ 52] = Gcz*cay1 + Ryay1;
        ERP[i][ 53] = 0.5*(caxxw - cayyw);
        ERP[i][ 54] = SQR3I*(cazzw - 0.5*(caxxw + cayyw));
        ERP[i][ 55] = Gcx*cby1         + Rczy1;
        ERP[i][ 56] = Gcz*cbx1 + Rybx1 + Rcxx1;
        ERP[i][ 57] = Gcz*cby1 + Ryby1 + Rcxy1;
        ERP[i][ 58] = 0.5*(cbxxw - cbyyw);
        ERP[i][ 59] = SQR3I*(cbzzw - 0.5*(cbxxw + cbyyw));
        ERP[i][ 60] = Gcx*ccy1;
        ERP[i][ 61] = Gcz*ccx1 + Rycx1 + Rcyx1;
        ERP[i][ 62] = Gcz*ccy1 + Rycy1 + Rcyy1;
        ERP[i][ 63] = 0.5*(ccxxw - ccyyw);
        ERP[i][ 64] = SQR3I*(cczzw - 0.5*(ccxxw + ccyyw));
        ERP[i][ 65] = Gcx*cdy1         + Rcxy1;
        ERP[i][ 66] = Gcz*cdx1 + Rydx1;
        ERP[i][ 67] = Gcz*cdy1 + Rydy1;
        ERP[i][ 68] = 0.5*(cdxxw - cdyyw);
        ERP[i][ 69] = SQR3I*(cdzzw - 0.5*(cdxxw + cdyyw));
        ERP[i][ 70] = Gcx*cey1         - SQR3I *Rcxy1;
        ERP[i][ 71] = Gcz*cex1 + Ryex1 + SQR3I2*Rczx1;
        ERP[i][ 72] = Gcz*cey1 + Ryey1 + SQR3I2*Rczy1;
        ERP[i][ 73] = 0.5*(cexxw - ceyyw);
        ERP[i][ 74] = SQR3I*(cezzw - 0.5*(cexxw + ceyyw));

        ERP[i][ 75] = Gcx*day1 + Rxay1 + Rdyy1;
        ERP[i][ 76] = Gcz*dax1;
        ERP[i][ 77] = Gcz*day1;
        ERP[i][ 78] = 0.5*(daxxw - dayyw);
        ERP[i][ 79] = SQR3I*(dazzw - 0.5*(daxxw + dayyw));
        ERP[i][ 80] = Gcx*dby1 + Rxby1 + Rdzy1;
        ERP[i][ 81] = Gcz*dbx1         + Rdxx1;
        ERP[i][ 82] = Gcz*dby1         + Rdxy1;
        ERP[i][ 83] = 0.5*(dbxxw - dbyyw);
        ERP[i][ 84] = SQR3I*(dbzzw - 0.5*(dbxxw + dbyyw));
        ERP[i][ 85] = Gcx*dcy1 + Rxcy1;
        ERP[i][ 86] = Gcz*dcx1         + Rdyx1;
        ERP[i][ 87] = Gcz*dcy1         + Rdyy1;
        ERP[i][ 88] = 0.5*(dcxxw - dcyyw);
        ERP[i][ 89] = SQR3I*(dczzw - 0.5*(dcxxw + dcyyw));
        ERP[i][ 90] = Gcx*ddy1 + Rxdy1 + Rdxy1;
        ERP[i][ 91] = Gcz*ddx1;
        ERP[i][ 92] = Gcz*ddy1;
        ERP[i][ 93] = 0.5*(ddxxw - ddyyw);
        ERP[i][ 94] = SQR3I*(ddzzw - 0.5*(ddxxw + ddyyw));
        ERP[i][ 95] = Gcx*dey1 + Rxey1 - SQR3I *Rdxy1;
        ERP[i][ 96] = Gcz*dex1         + SQR3I2*Rdzx1;
        ERP[i][ 97] = Gcz*dey1         + SQR3I2*Rdzy1;
        ERP[i][ 98] = 0.5*(dexxw - deyyw);
        ERP[i][ 99] = SQR3I*(dezzw - 0.5*(dexxw + deyyw));

        ERP[i][100] = Gcx*eay1 - SQR3I *Rxay1 + Reyy1;
        ERP[i][101] = Gcz*eax1 + SQR3I2*Rzax1;
        ERP[i][102] = Gcz*eay1 + SQR3I2*Rzay1;
        ERP[i][103] = 0.5*(eaxxw - eayyw);
        ERP[i][104] = SQR3I*(eazzw - 0.5*(eaxxw + eayyw));
        ERP[i][105] = Gcx*eby1 - SQR3I *Rxby1 + Rezy1;
        ERP[i][106] = Gcz*ebx1 + SQR3I2*Rzbx1 + Rexx1;
        ERP[i][107] = Gcz*eby1 + SQR3I2*Rzby1 + Rexy1;
        ERP[i][108] = 0.5*(ebxxw - ebyyw);
        ERP[i][109] = SQR3I*(ebzzw - 0.5*(ebxxw + ebyyw));
        ERP[i][110] = Gcx*ecy1 - SQR3I *Rxcy1;
        ERP[i][111] = Gcz*ecx1 + SQR3I2*Rzcx1 + Reyx1;
        ERP[i][112] = Gcz*ecy1 + SQR3I2*Rzcy1 + Reyy1;
        ERP[i][113] = 0.5*(ecxxw - ecyyw);
        ERP[i][114] = SQR3I*(eczzw - 0.5*(ecxxw + ecyyw));
        ERP[i][115] = Gcx*edy1 - SQR3I *Rxdy1 + Rexy1;
        ERP[i][116] = Gcz*edx1 + SQR3I2*Rzdx1;
        ERP[i][117] = Gcz*edy1 + SQR3I2*Rzdy1;
        ERP[i][118] = 0.5*(edxxw - edyyw);
        ERP[i][119] = SQR3I*(edzzw - 0.5*(edxxw + edyyw));
        ERP[i][120] = Gcx*eey1 - SQR3I *(Rxey1 + Rexy1);
        ERP[i][121] = Gcz*eex1 + SQR3I2*(Rzex1 + Rezx1);
        ERP[i][122] = Gcz*eey1 + SQR3I2*(Rzey1 + Rezy1);
        ERP[i][123] = 0.5*(eexxw - eeyyw);
        ERP[i][124] = SQR3I*(eezzw - 0.5*(eexxw + eeyyw));
    }
}


bool DfEri::isCutoff_IK(const std::size_t orbA, const std::size_t orbC) const
{
    bool bAnswer = true;

    // distance
    //const int nuA = this->IshellTable[1][nIShell];
    //const int nuC = this->KshellTable[1][nKShell];
    //const double AC2 = this->coord_[nuA].squareDistanceFrom(this->coord_[nuC]);
    //const double AC = std::sqrt(AC2);

    // zeta
    //const int iAOs = IshellTable[5][nIShell];
    //const std::size_t orbIndexA = this->IDTableO[iAOs];
    const double dZetaA = this->pOrbitalInfo_->getMinExponent(orbA);
    const TlPosition A = this->pOrbitalInfo_->getPosition(orbA);

    //const std::size_t kAOs = KshellTable[4][nKShell];
    //const std::size_t orbC = this->IDTableD[kAOs];
    const double dZetaC = this->pOrbitalInfo_Density_->getExponent(orbC, 0);
    const TlPosition C = this->pOrbitalInfo_Density_->getPosition(orbC);

    const double AC2 = A.squareDistanceFrom(C);
    const double AC = std::sqrt(AC2);

    // R
    const double dZeta = dZetaA + dZetaC;
    const double dZetaAC = dZetaA * dZetaC;
    const double t = dZetaAC * AC2 / dZeta;

    TlFmt& FmT = TlFmt::getInstance();
    const double erfc = FmT.getFmT0(0, t);
    const double inv_erfc = 1.0 / erfc;

    const double R1 = std::sqrt(2.0 * dZetaA) * inv_erfc +1.0;
    const double R2 = std::sqrt(2.0 * dZetaC) * inv_erfc +1.0;

    double dJudge = (R1 + R2 +2.0) - AC;

    if (dJudge > this->m_dCutValue) {
        bAnswer = false;
    }

    return bAnswer;
}

