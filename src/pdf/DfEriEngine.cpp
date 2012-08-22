#include <algorithm>
#include <iostream>
#include "DfEriEngine.h"
#include "TlFmt.h"
#include "TlUtils.h"
#include "TlMath.h"

// デバッグ用フラグ
// #define DEBUG_HGP
// #define DEBUG_CALC_R0
// #define DEBUG_CHOICE
// #define DEBUG_R_DASH
// #define DEBUG_R_2DASH
// #define DEBUG_CONTRACT_BRA
// #define DEBUG_BRA_ERI
// #define DEBUG_CALC_PQ
// #define DEBUG_EQ43
// #define DEBUG_EQ44
// #define DEBUG_EQ45
// #define DEBUG_EQ46
// #define DEBUG_EQ47
// #define DEBUG_OUTPUT
// #define DEBUG_TRANSFORM_6D_TO_5D

#define USE_CACHED_ROUTE // choice()をキャッシングする

/// 結果出力用のサイズ
/// d軌道の6Dから5Dへの変換領域にも利用するため、
/// d軌道のサイズを6Dベースとして確保している。
/// = 6(d) * 6 * 6 * 6 * 3
const int DfEriEngine::OUTPUT_BUFFER_SIZE = 3 * 6 * 6 * 6 * 6;

const double DfEriEngine::INV_SQRT3 = 1.0 / sqrt(3.0);

const TlAngularMomentumVector DfEriEngine::E1_[3] = {
    TlAngularMomentumVector(1, 0, 0),
    TlAngularMomentumVector(0, 1, 0),
    TlAngularMomentumVector(0, 0, 1)
};

const TlAngularMomentumVector DfEriEngine::E2_[3] = {
    TlAngularMomentumVector(2, 0, 0),
    TlAngularMomentumVector(0, 2, 0),
    TlAngularMomentumVector(0, 0, 2)
};


////////////////////////////////////////////////////////////////////////////////
DfEriEngine::DfEriEngine()
    : log_(TlLogging::getInstance())
{
    this->p0M_ = new double*[ERI_L_MAX];
    for (int i = 0; i < ERI_L_MAX; ++i) {
        this->p0M_[i] = new double[ERI_KPKQ_MAX];
    }

    this->pRM_ = new double*[ERI_NUM_OF_RM_KINDS];
    for (int i = 0; i < ERI_NUM_OF_RM_KINDS; ++i) {
        this->pRM_[i] = new double[ERI_KPKQ_MAX];
    }

    this->pContractBraCoef_ = new double[ERI_KPKQ_MAX];
    
    // abp(r)cdq
    const std::size_t P_PRIME_MAX3 = ERI_P_PRIME_MAX * ERI_P_PRIME_MAX * ERI_P_PRIME_MAX;
    const std::size_t max_abpRcdq_index = P_PRIME_MAX3 * P_PRIME_MAX3 * ERI_NUM_OF_R_KINDS;
    this->p_abpRcdq_ = new double[max_abpRcdq_index];

    // ERIメモリ確保
    this->ERI_bra_.resize(ERI_NUM_OF_ERI_STATES);
    this->ERI_ket_.resize(ERI_NUM_OF_ERI_STATES);
    // for (int i = 0; i < ERI_L_MAX; ++i) {
    //     const TlAngularMomentumVectorSet amvs(i);
    //     amvs.size();
    // }

    
    this->pGmT_ = new double[ERI_L_MAX +1];
    this->pE4_ = new E4[ERI_KPKQ_MAX];

    this->nR_dash_.resize(ERI_NR_DASH_SIZE);
    
    this->WORK = new double[OUTPUT_BUFFER_SIZE];
    this->WORK_A = new double[OUTPUT_BUFFER_SIZE];
    this->WORK_B = new double[OUTPUT_BUFFER_SIZE];
    this->WORK_C = new double[OUTPUT_BUFFER_SIZE];
    this->WORK_D = new double[OUTPUT_BUFFER_SIZE];
    
    // cutoff
    this->primitiveLevelThreshold_ = 1.0E-20;

#ifdef CHECK_MAX_COUNT
    this->maxSizeOf_sumOfAngularMomentums_ = 0;
    this->max_a_bar_ = 0;
    this->max_b_bar_ = 0;
    this->max_a_ = 0;
    this->max_b_ = 0;
    this->max_p_ = 0;
    this->max_a_prime_ = 0;
    this->max_b_prime_ = 0;
    this->max_p_prime_ = 0;
    this->maxSizeOf_nR_dash_ = 0;
    this->maxNumOfAMVs_ = 0;
    this->maxERI_batch_ = 0;
#endif // CHECK_MAX_COUNT
}


DfEriEngine::~DfEriEngine()
{
    if (this->p0M_ != NULL) {
        for (int i = 0; i < ERI_L_MAX; ++i) {
            if (this->p0M_[i] != NULL) {
                delete[] this->p0M_[i];
                this->p0M_[i] = NULL;
            }
        }
        delete[] this->p0M_;
        this->p0M_ = NULL;
    }

    if (this->pRM_ != NULL) {
        for (int i = 0; i < ERI_NUM_OF_RM_KINDS; ++i) {
            if (this->pRM_[i] != NULL) {
                delete[] this->pRM_[i];
                this->pRM_[i] = NULL;
            }
        }
        delete[] this->pRM_;
        this->pRM_ = NULL;
    }

    if (this->pContractBraCoef_ != NULL) {
        delete[] this->pContractBraCoef_;
        this->pContractBraCoef_ = NULL;
    }

    if (this->p_abpRcdq_ != NULL) {
        delete[] this->p_abpRcdq_;
        this->p_abpRcdq_ = NULL;
    }

    delete[] this->pGmT_;
    this->pGmT_ = NULL;

    delete[] this->pE4_;
    this->pE4_ = NULL;
    
    delete[] this->WORK;
    this->WORK = NULL;
    delete[] this->WORK_A;
    this->WORK_A = NULL;
    delete[] this->WORK_B;
    this->WORK_B = NULL;
    delete[] this->WORK_C;
    this->WORK_C = NULL;
    delete[] this->WORK_D;
    this->WORK_D = NULL;

#ifdef CHECK_MAX_COUNT
    std::cerr << ">>>> DfEriEngine constants" << std::endl;
    std::cerr << TlUtils::format("max sumOfAngularMomentums = %ld", this->maxSizeOf_sumOfAngularMomentums_)
              << std::endl;
    std::cerr << TlUtils::format("max (a~, b~, a, b, p, a', b', p') = (%d, %d, %d, %d, %d, %d, %d, %d)",
                                 this->max_a_bar_, this->max_b_bar_,
                                 this->max_a_, this->max_b_, this->max_p_,
                                 this->max_a_prime_, this->max_b_prime_, this->max_p_prime_)
              << std::endl;
    std::cerr << TlUtils::format("max nR_dash_ = %ld", this->maxSizeOf_nR_dash_)
              << std::endl;
    std::cerr << TlUtils::format("max num of AMVs = %d", this->maxNumOfAMVs_)
              << std::endl;
    std::cerr << TlUtils::format("max ERI batch = %d", this->maxERI_batch_)
              << std::endl;
    std::cerr << std::endl;
#endif // CHECK_MAX_COUNT
}


DfEriEngine::CGTO_Pair DfEriEngine::getCGTO_pair(const TlOrbitalInfoObject& orbInfo,
                                                 const index_type shellIndexP,
                                                 const index_type shellIndexQ,
                                                 const double cutoffThreshold)
{
    const TlPosition P = orbInfo.getPosition(shellIndexP);
    const int shellTypeP = orbInfo.getShellType(shellIndexP);
    int numOfContractionsP = orbInfo.getCgtoContraction(shellIndexP);

    TlPosition Q(0.0, 0.0, 0.0);
    int shellTypeQ = 0;
    int numOfContractionsQ = 1;
    if (shellIndexQ >= 0) {
        Q = orbInfo.getPosition(shellIndexQ);
        shellTypeQ = orbInfo.getShellType(shellIndexQ);
        numOfContractionsQ = orbInfo.getCgtoContraction(shellIndexQ);
    }
    const TlPosition PQ = Q - P;

    // cutoff
    // TlOrbitalInfoObjectオブジェクトは指数の大きい順に並んでいるので
    // Pの1番目(最大)とQを比較すれば良い。
    // static const double INV_EQ32_COEF = 1.0 / (std::pow(2.0 * TlMath::PI(), 0.25) * TlMath::PI());
    // const double threshold = cutoffThreshold * INV_EQ32_COEF;
    // const double distance2 = PQ.squareDistanceFrom();
    // if (shellIndexQ >= 0) {
    //     const double exponentQ = orbInfo.getExponent(shellIndexQ, 0);
    //     for (int i = 0; i < numOfContractionsP; ++i) {
    //         const double exponentP = orbInfo.getExponent(shellIndexP, i);
    //         const double zeta = exponentP + exponentQ;
    //         const double inv_zeta = 1.0 / zeta;
    //         const double Kab = std::fabs(inv_zeta * std::exp(- exponentP * exponentQ * inv_zeta * distance2));
    //         if (Kab < threshold) {
    //             numOfContractionsP = i +1;
    //             break;
    //         }
    //     }

    //     const double exponentP = orbInfo.getExponent(shellIndexP, 0);
    //     for (int i = 0; i < numOfContractionsQ; ++i) {
    //         const double exponentQ = orbInfo.getExponent(shellIndexQ, i);
    //         const double zeta = exponentP + exponentQ;
    //         const double inv_zeta = 1.0 / zeta;
    //         const double Kab = std::fabs(inv_zeta * std::exp(- exponentP * exponentQ * inv_zeta * distance2));
    //         if (Kab < threshold) {
    //             numOfContractionsQ = i +1;
    //             break;
    //         }
    //     }
    // }

    // make pair
    DfEriEngine::PGTO_Pairs pgtosPQ(numOfContractionsP * numOfContractionsQ);
    int index = 0;
    for (int pP = 0; pP < numOfContractionsP; ++pP) {
        const double coefP = orbInfo.getCoefficient(shellIndexP, pP);
        const double expP = orbInfo.getExponent(shellIndexP, pP);
        
        for (int pQ = 0; pQ < numOfContractionsQ; ++pQ) {
            double coefQ = 1.0;
            double expQ  = 0.0;
            if (shellIndexQ >= 0) {
                coefQ = orbInfo.getCoefficient(shellIndexQ, pQ);
                expQ = orbInfo.getExponent(shellIndexQ, pQ);
            }
            
            DfEriEngine::PGTO_Pair pgtoPair(shellTypeP, coefP, expP,
                                            shellTypeQ, coefQ, expQ,
                                            P, Q);
            pgtosPQ[index] = pgtoPair;
            ++index;
        }
    }
    CGTO_Pair pairPQ(PQ, pgtosPQ);
    std::sort(pairPQ.PS.begin(), pairPQ.PS.end(), PGTO_sort_functor_cmp());

    return pairPQ;
}


void DfEriEngine::calc(const Query& qAB, const Query& qCD,
                       const CGTO_Pair& IJ, const CGTO_Pair& KL)
{
    assert((0 <= qAB.a_bar) && (qAB.a_bar < ERI_A_BAR_MAX));
    assert((0 <= qAB.b_bar) && (qAB.b_bar < ERI_B_BAR_MAX));
    assert((0 <= qAB.a) && (qAB.a < ERI_A_MAX));
    assert((0 <= qAB.b) && (qAB.b < ERI_B_MAX));
    
    // this->time_calc_all_.start();
    // initialize
    this->sumOfAngularMomentums_ = qAB.sum() + qCD.sum();
#ifdef CHECK_MAX_COUNT
    this->maxSizeOf_sumOfAngularMomentums_ = std::max(this->maxSizeOf_sumOfAngularMomentums_,
                                                      this->sumOfAngularMomentums_);
#endif // CHECK_MAX_COUNT

    this->initialize();
    // fill is too slow!
    std::fill_n(this->WORK, OUTPUT_BUFFER_SIZE, 0.0); 
    // for (int i = 0; i < OUTPUT_BUFFER_SIZE; ++i) {
    //     this->WORK[i] = 0.0;
    // }
    
    // 総角運動量を求める
    // const int a = qAB.a;
    // const int b = qAB.b;
    // const int a_bar = qAB.a_bar;
    // const int b_bar = qAB.b_bar;
    // const int c = qCD.a;
    // const int d = qCD.b;
    // const int c_bar = qCD.a_bar;
    // const int d_bar = qCD.b_bar;
    // std::cerr << TlUtils::format("calc (%d %d, %d %d|%d %d, %d %d)=%d",
    //                              a_bar, b_bar, a, b,
    //                              c_bar, d_bar, c, d, this->L_) << std::endl;
    const int numOfPGTOs_IJ = IJ.PS.size();
    const int numOfPGTOs_KL = KL.PS.size();

    this->AB_ = IJ.AB;
    this->CD_ = KL.AB;

#ifdef DEBUG_HGP
    // for debug
    {
        const std::string AB_str = TlUtils::format("[%f %f %f]", this->AB_[0], this->AB_[1], this->AB_[2]);
        const std::string CD_str = TlUtils::format("[%f %f %f]", this->CD_[0], this->CD_[1], this->CD_[2]);
        std::cerr << TlUtils::format(">>>> BEGIN\n(%d %d, %d %d |%d %d, %d %d) AB=%s, CD=%s",
                                     qAB.a_bar, qAB.b_bar, qAB.a, qAB.b,
                                     qCD.a_bar, qCD.b_bar, qCD.a, qCD.b,
                                     AB_str.c_str(), CD_str.c_str())
                  << std::endl;
    }
#endif // DEBUG_HGP

    
    // 一度に計算するbra-ketに切り分ける
    this->bra_ = IJ.PS;
    
    PGTO_Pairs ket;
    const int block_KL = (ERI_KPKQ_MAX -1) / numOfPGTOs_IJ;
    const div_t blocks = std::div(numOfPGTOs_KL, block_KL);
    ket.resize(block_KL);
    for (int block = 0; block < blocks.quot; ++block) {
        const int base = block_KL * block;
        std::copy(KL.PS.begin() +base, KL.PS.begin() +base +block_KL, ket.begin());
        this->ket_ = ket;

#ifdef DEBUG_HGP
        // for debug
        std::cerr << TlUtils::format("block1: KP=%d, KQ=%d",
                                     this->bra_.size(), this->ket_.size())
                  << std::endl;
#endif // DEBUG_HGP
        
        this->calc(qAB, qCD);
        this->copyResultsToOutputBuffer(qAB, qCD, this->WORK);
    }

    {
        const int base = block_KL * blocks.quot;
        ket.resize(blocks.rem);
        std::copy(KL.PS.begin() +base, KL.PS.begin() +base +blocks.rem, ket.begin());
        this->ket_ = ket;

#ifdef DEBUG_HGP
        // for debug
        std::cerr << TlUtils::format("block2: KP=%d, KQ=%d",
                                     this->bra_.size(), this->ket_.size())
                  << std::endl;
#endif // DEBUG_HGP
        
        this->calc(qAB, qCD);
        this->copyResultsToOutputBuffer(qAB, qCD, this->WORK);
    }

    //for debug
    // {
    //     std::cerr << TlUtils::format(">>>> 6D >>>> (%d %d %d %d | %d %d %d %d)",
    //                                  a_bar, b_bar, c_bar, d_bar, a, b, c, d)
    //               << std::endl;
    //     int A_MAX = a_bar * (a_bar + 3) / 2 + 1;
    //     int B_MAX = b_bar * (b_bar + 3) / 2 + 1;
    //     int C_MAX = c_bar * (c_bar + 3) / 2 + 1;
    //     int D_MAX = d_bar * (d_bar + 3) / 2 + 1;
    //     int AMAX = a * (a + 3) / 2 + 1;
    //     int BMAX = b * (b + 3) / 2 + 1;
    //     int CMAX = c * (c + 3) / 2 + 1;
    //     int DMAX = d * (d + 3) / 2 + 1;
    //     int index = 0;
    //     for (int A_ = 0; A_ < A_MAX; ++A_) {
    //         for (int B_ = 0; B_ < B_MAX; ++B_) {
    //             for (int C_ = 0; C_ < C_MAX; ++C_) {
    //                 for (int D_ = 0; D_ < D_MAX; ++D_) {
    //                     for (int A = 0; A < AMAX; ++A) {
    //                         for (int B = 0; B < BMAX; ++B) {
    //                             for (int C = 0; C < CMAX; ++C) {
    //                                 for (int D = 0; D < DMAX; ++D) {
    //                                     double value = this->WORK[index];
    //                                     if (std::fabs(value) > 1.0E-20) {
    //                                         std::cerr << TlUtils::format("CART %4d (%d %d, %d %d | %d %d, %d %d) = % f",
    //                                                                      index,
    //                                                                      A_, B_, A, B,
    //                                                                      C_, D_, C, D,
    //                                                                      value)
    //                                                   << std::endl;
    //                                     }
    //                                     ++index;
    //                                 }
    //                             }
    //                         }
    //                     }
    //                 }
    //             }
    //         }
    //     }
    //     std::cerr << std::endl;
    // }
    
    // transform 6D to 5D
    this->transform6Dto5D(qAB, qCD, this->WORK);
    
   
    // this->time_calc_all_.stop();
#ifdef DEBUG_HGP
    std::cerr << "<<<<END\n" << std::endl;
#endif // DEBUG_HGP
}


void DfEriEngine::calcGrad(const Query& in_qAB, const Query& in_qCD,
                           const CGTO_Pair& IJ, const CGTO_Pair& KL)
{
    // initialize
    this->sumOfAngularMomentums_ = in_qAB.a + in_qAB.b + in_qCD.a + in_qCD.b + 1;
#ifdef CHECK_MAX_COUNT
    this->maxSizeOf_sumOfAngularMomentums_ = std::max(this->maxSizeOf_sumOfAngularMomentums_,
                                                      this->sumOfAngularMomentums_);
#endif // CHECK_MAX_COUNT
    
    this->initialize();
    std::fill(this->WORK_A, this->WORK_A + OUTPUT_BUFFER_SIZE, 0.0); 
    std::fill(this->WORK_B, this->WORK_B + OUTPUT_BUFFER_SIZE, 0.0); 
    std::fill(this->WORK_C, this->WORK_C + OUTPUT_BUFFER_SIZE, 0.0); 
    
    const int numOfPGTOs_IJ = IJ.PS.size();
    const int numOfPGTOs_KL = KL.PS.size();
    this->AB_ = IJ.AB;
    this->CD_ = KL.AB;

    {
        const Query qAB00(0, 0, in_qAB.a, in_qAB.b);
        const Query qCD00(0, 0, in_qCD.a, in_qCD.b);

        // 一度に計算するbra-ketに切り分ける
        this->bra_ = IJ.PS;
        
        PGTO_Pairs ket;
        const int block_KL = (ERI_KPKQ_MAX -1) / numOfPGTOs_IJ;
        const div_t blocks = std::div(numOfPGTOs_KL, block_KL);
        ket.resize(block_KL);
        for (int block = 0; block < blocks.quot; ++block) {
            const int base = block_KL * block;
            std::copy(KL.PS.begin() +base, KL.PS.begin() +base +block_KL, ket.begin());
            this->ket_ = ket;
            
            this->calcGrad(qAB00, qCD00);
            //this->copyResultsToOutputBuffer(qAB10, qCD00, this->WORK_A);
        }
        {
            const int base = block_KL * blocks.quot;
            ket.resize(blocks.rem);
            std::copy(KL.PS.begin() +base, KL.PS.begin() +base +blocks.rem, ket.begin());
            this->ket_ = ket;
            
            this->calcGrad(qAB00, qCD00);
            //this->copyResultsToOutputBuffer(qAB10, qCD00, this->WORK_A);
        }

        // transform 6D to 5D
        const Query qAB10(1, 0, in_qAB.a, in_qAB.b);
        const Query qAB01(0, 1, in_qAB.a, in_qAB.b);
        const Query qCD10(1, 0, in_qCD.a, in_qCD.b);
        this->transform6Dto5D(qAB10, qCD00, this->WORK_A);
        this->transform6Dto5D(qAB01, qCD00, this->WORK_B);
        this->transform6Dto5D(qAB00, qCD10, this->WORK_C);
    }
    
    // D'
    {
        const Query qAB00(0, 0, in_qAB.a, in_qAB.b);
        const Query qCD01(0, 1, in_qCD.a, in_qCD.b);
        this->compD(qAB00, qCD01);
    }
}


void DfEriEngine::initialize()
{
    const int sumOfAngularMomentums = this->sumOfAngularMomentums_;
    const int KPKQ = this->bra_.size() * this->ket_.size();

    // p0M_
    for (int m = 0; m <= sumOfAngularMomentums; ++m) {
        std::fill(this->p0M_[m], this->p0M_[m] + KPKQ, (double)0.0);
    }

    // pRM_
    for (int r = 0; r <= sumOfAngularMomentums; ++r) {
        const int indexR0 = this->indexRM(r, 0);
        std::fill(this->pRM_[indexR0], this->pRM_[indexR0] + KPKQ, 0.0);
    }
}


void DfEriEngine::copyResultsToOutputBuffer(const DfEriEngine::Query& qAB,
                                            const DfEriEngine::Query& qCD,
                                            double* pOutput)
{
    const int a_bar = qAB.a_bar;
    const int b_bar = qAB.b_bar;
    const int a = qAB.a;
    const int b = qAB.b;
    const int c_bar = qCD.a_bar;
    const int d_bar = qCD.b_bar;
    const int c = qCD.a;
    const int d = qCD.b;

    const ERI_State ket_state(c_bar, d_bar, c, d, 0, 0, 0, 0);
    const std::size_t ketStateIndex = ket_state.index();
    
    const TlAngularMomentumVectorSet amvsAbar(a_bar);
    const TlAngularMomentumVectorSet amvsBbar(b_bar);
    const TlAngularMomentumVectorSet amvsCbar(c_bar);
    const TlAngularMomentumVectorSet amvsDbar(d_bar);
    const TlAngularMomentumVectorSet amvsA(a);
    const TlAngularMomentumVectorSet amvsB(b);
    const TlAngularMomentumVectorSet amvsC(c);
    const TlAngularMomentumVectorSet amvsD(d);
    const TlAngularMomentumVector amvQ(0, 0, 0);

    const int numOfAmvsAbar = amvsAbar.size();
    const int numOfAmvsBbar = amvsBbar.size();
    const int numOfAmvsCbar = amvsCbar.size();
    const int numOfAmvsDbar = amvsDbar.size();
    const int numOfAmvsA = amvsA.size();
    const int numOfAmvsB = amvsB.size();
    const int numOfAmvsC = amvsC.size();
    const int numOfAmvsD = amvsD.size();
    
    int index = 0;
    for (int iAbar = 0; iAbar < numOfAmvsAbar; ++iAbar) {
        const TlAngularMomentumVector amvAbar = amvsAbar.get(iAbar);

        for (int iBbar = 0; iBbar < numOfAmvsBbar; ++iBbar) {
            const TlAngularMomentumVector amvBbar = amvsBbar.get(iBbar);

            for (int iCbar = 0; iCbar < numOfAmvsCbar; ++iCbar) {
                const TlAngularMomentumVector amvCbar = amvsCbar.get(iCbar);
                
                for (int iDbar = 0; iDbar < numOfAmvsDbar; ++iDbar) {
                    const TlAngularMomentumVector amvDbar = amvsDbar.get(iDbar);

                    for (int iA = 0; iA < numOfAmvsA; ++iA) {
                        const TlAngularMomentumVector amvA = amvsA.get(iA);
                
                        for (int iB = 0; iB < numOfAmvsB; ++iB) {
                            const TlAngularMomentumVector amvB = amvsB.get(iB);
                            const int bra_index = ((iAbar*numOfAmvsBbar +iBbar)*numOfAmvsA +iA)*numOfAmvsB +iB;
                    
                            for (int iC = 0; iC < numOfAmvsC; ++iC) {
                                const TlAngularMomentumVector amvC = amvsC.get(iC);

                                for (int iD = 0; iD < numOfAmvsD; ++iD) {
                                    const TlAngularMomentumVector amvD = amvsD.get(iD);
                                    const int ket_index = this->index(amvCbar, amvDbar,
                                                                      amvC, amvD, amvQ);
                                    pOutput[index] += this->ERI_ket_[ketStateIndex][ket_index][bra_index];
                                    ++index;

#ifdef DEBUG_OUTPUT
                                    // for debug
                                    {
                                        double v = eri[ket_state][ket_index][bra_index];
                                        std::cerr << TlUtils::format("output[%3d] (%s %s, %s %s|%s %s, %s %s) = %f",
                                                                     index,
                                                                     amv_a_bar.debugOut().c_str(),
                                                                     amv_b_bar.debugOut().c_str(),
                                                                     amv_a.debugOut().c_str(),
                                                                     amv_b.debugOut().c_str(),
                                                                     amv_c_bar.debugOut().c_str(),
                                                                     amv_d_bar.debugOut().c_str(),
                                                                     amv_c.debugOut().c_str(),
                                                                     amv_d.debugOut().c_str(),
                                                                     v)
                                                  << std::endl;
                                    }
#endif //DEBUG_OUTPUT
                                } // d
                            } // c
                        } // b
                    } // a
                }
            }
        }
    }
}


void DfEriEngine::transform6Dto5D(const DfEriEngine::Query& qAB,
                                  const DfEriEngine::Query& qCD,
                                  double* pOutput)
{
    const int a_bar = qAB.a_bar;
    const int b_bar = qAB.b_bar;
    const int c_bar = qCD.a_bar;
    const int d_bar = qCD.b_bar;
    const int a = qAB.a;
    const int b = qAB.b;
    const int c = qCD.a;
    const int d = qCD.b;
    
    if ((a < 2) && (b < 2) && (c < 2) && (d < 2)) {
        // nothing to do
        return;
    }

    double pBuf[OUTPUT_BUFFER_SIZE];
    
    // 6D ベースの要素数(1, 3, 6, 10, ...)
    // 5D ベースの場合は(2n +1: 1, 3, 5, 7, ...)
    int I_ = a_bar * (a_bar + 3) / 2 + 1;
    int J_ = b_bar * (b_bar + 3) / 2 + 1;
    int K_ = c_bar * (c_bar + 3) / 2 + 1;
    int L_ = d_bar * (d_bar + 3) / 2 + 1;
    int I = a * (a + 3) / 2 + 1;
    int J = b * (b + 3) / 2 + 1;
    int K = c * (c + 3) / 2 + 1;
    int L = d * (d + 3) / 2 + 1;

#ifdef DEBUG_TRANSFORM_6D_TO_5D
    std::cerr << TlUtils::format("DfEriEngine::transform6Dto5D() in: a=%d, b=%d, c=%d, d=%d",
                                 a, b, c, d)
              << std::endl;
#endif // DEBUG_TRANSFORM_6D_TO_5D
    
    if (a == 2) {
        this->transform6Dto5D_i(I_, J_, K_, L_, J, K, L, pOutput, pBuf);
        I = 5;
        const std::size_t end = I_ * J_ * K_ * L_ * I * J * K * L;
        std::copy(pBuf, pBuf + end, pOutput);
    }
    if (b == 2) {
        this->transform6Dto5D_j(I_, J_, K_, L_, I, K, L, pOutput, pBuf);
        J = 5;
        const std::size_t end = I_ * J_ * K_ * L_ * I * J * K * L;
        std::copy(pBuf, pBuf + end, pOutput);
    }
    if (c == 2) {
        this->transform6Dto5D_k(I_, J_, K_, L_, I, J, L, pOutput, pBuf);
        K = 5;
        const std::size_t end = I_ * J_ * K_ * L_ * I * J * K * L;
        std::copy(pBuf, pBuf + end, pOutput);
    }
    if (d == 2) {
        this->transform6Dto5D_l(I_, J_, K_, L_, I, J, K, pOutput, pBuf);
        L = 5;
        const std::size_t end = I_ * J_ * K_ * L_ * I * J * K * L;
        std::copy(pBuf, pBuf + end, pOutput);
    }
}


void DfEriEngine::transform6Dto5D_i(const int I_, const int J_, const int K_, const int L_,
                              const int J, const int K, const int L,
                              const double* pInput, double* pOutput)
{
#ifdef DEBUG_TRANSFORM_6D_TO_5D
    std::cerr << TlUtils::format("DfEriEngine::transform6Dto5D_i() in: I^=%d, J^=%d, K^=%d, L^=%d, J=%d, K=%d, L=%d",
                                 I_, J_, K_, L_, J, K, L)
              << std::endl;
#endif // DEBUG_TRANSFORM_6D_TO_5D
    
    for (int i_ = 0; i_ < I_; ++i_) {
        for (int j_ = 0; j_ < J_; ++j_) {
            for (int k_ = 0; k_ < K_; ++k_) {
                for (int l_ = 0; l_ < L_; ++l_) {

                    for (int j = 0; j < J; ++j) {
                        for (int k = 0; k < K; ++k) {
                            for (int l = 0; l < L; ++l) {
                                const double xx = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ + l_)*6 +0)*J +j)*K +k)*L +l];
                                const double xy = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ + l_)*6 +1)*J +j)*K +k)*L +l];
                                const double xz = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ + l_)*6 +2)*J +j)*K +k)*L +l];
                                const double yy = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ + l_)*6 +3)*J +j)*K +k)*L +l];
                                const double yz = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ + l_)*6 +4)*J +j)*K +k)*L +l];
                                const double zz = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ + l_)*6 +5)*J +j)*K +k)*L +l];

                                const int xy_5d   = ((((((i_*J_ +j_)*K_ +k_)*L_ + l_)*5 +0)*J +j)*K +k)*L +l;
                                const int xz_5d   = ((((((i_*J_ +j_)*K_ +k_)*L_ + l_)*5 +1)*J +j)*K +k)*L +l;
                                const int yz_5d   = ((((((i_*J_ +j_)*K_ +k_)*L_ + l_)*5 +2)*J +j)*K +k)*L +l;
                                const int xxyy_5d = ((((((i_*J_ +j_)*K_ +k_)*L_ + l_)*5 +3)*J +j)*K +k)*L +l;
                                const int rr_5d   = ((((((i_*J_ +j_)*K_ +k_)*L_ + l_)*5 +4)*J +j)*K +k)*L +l;
                                pOutput[xy_5d] = xy;
                                pOutput[xz_5d] = xz;
                                pOutput[yz_5d] = yz;
                                pOutput[xxyy_5d] = 0.5*(xx-yy);
                                pOutput[rr_5d] = INV_SQRT3 * (zz - 0.5 * (xx + yy));
                            }
                        }
                    }
                }
            }
        }
    }
}


void DfEriEngine::transform6Dto5D_j(const int I_, const int J_, const int K_, const int L_,
                              const int I, const int K, const int L,
                              const double* pInput, double* pOutput)
{
#ifdef DEBUG_TRANSFORM_6D_TO_5D
    std::cerr << TlUtils::format("DfEriEngine::transform6Dto5D_j() in: I^=%d, J^=%d, K^=%d, L^=%d, I=%d, K=%d, L=%d",
                                 I_, J_, K_, L_, I, K, L)
              << std::endl;
#endif // DEBUG_TRANSFORM_6D_TO_5D
    
    for (int i_ = 0; i_ < I_; ++i_) {
        for (int j_ = 0; j_ < J_; ++j_) {
            for (int k_ = 0; k_ < K_; ++k_) {
                for (int l_ = 0; l_ < L_; ++l_) {

                    for (int i = 0; i < I; ++i) {
                        for (int k = 0; k < K; ++k) {
                            for (int l = 0; l < L; ++l) {

                                const double xx = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ + l_)*I +i)*6 +0)*K +k)*L +l];
                                const double xy = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ + l_)*I +i)*6 +1)*K +k)*L +l];
                                const double xz = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ + l_)*I +i)*6 +2)*K +k)*L +l];
                                const double yy = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ + l_)*I +i)*6 +3)*K +k)*L +l];
                                const double yz = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ + l_)*I +i)*6 +4)*K +k)*L +l];
                                const double zz = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ + l_)*I +i)*6 +5)*K +k)*L +l];

                                const int xy_5d   = ((((((i_*J_ +j_)*K_ +k_)*L_ + l_)*I +i)*5 +0)*K +k)*L +l;
                                const int xz_5d   = ((((((i_*J_ +j_)*K_ +k_)*L_ + l_)*I +i)*5 +1)*K +k)*L +l;
                                const int yz_5d   = ((((((i_*J_ +j_)*K_ +k_)*L_ + l_)*I +i)*5 +2)*K +k)*L +l;
                                const int xxyy_5d = ((((((i_*J_ +j_)*K_ +k_)*L_ + l_)*I +i)*5 +3)*K +k)*L +l;
                                const int rr_5d   = ((((((i_*J_ +j_)*K_ +k_)*L_ + l_)*I +i)*5 +4)*K +k)*L +l;
                                pOutput[xy_5d] = xy;
                                pOutput[xz_5d] = xz;
                                pOutput[yz_5d] = yz;
                                pOutput[xxyy_5d] = 0.5*(xx-yy);
                                pOutput[rr_5d] = INV_SQRT3 * (zz - 0.5 * (xx + yy));
                            }
                        }
                    }
                }
            }
        }
    }
}


void DfEriEngine::transform6Dto5D_k(const int I_, const int J_, const int K_, const int L_,
                              const int I, const int J, const int L,
                              const double* pInput, double* pOutput)
{
#ifdef DEBUG_TRANSFORM_6D_TO_5D
    std::cerr << TlUtils::format("DfEriEngine::transform6Dto5D_i() in: I^=%d, J^=%d, K^=%d, L^=%d, I=%d, J=%d, L=%d",
                                 I_, J_, K_, L_, I, J, L)
              << std::endl;
#endif // DEBUG_TRANSFORM_6D_TO_5D
    
    for (int i_ = 0; i_ < I_; ++i_) {
        for (int j_ = 0; j_ < J_; ++j_) {
            for (int k_ = 0; k_ < K_; ++k_) {
                for (int l_ = 0; l_ < L_; ++l_) {

                    for (int i = 0; i < I; ++i) {
                        for (int j = 0; j < J; ++j) {
                            for (int l = 0; l < L; ++l) {
                                
                                const double xx = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ + l_)*I +i)*J +j)*6 +0)*L +l];
                                const double xy = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ + l_)*I +i)*J +j)*6 +1)*L +l];
                                const double xz = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ + l_)*I +i)*J +j)*6 +2)*L +l];
                                const double yy = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ + l_)*I +i)*J +j)*6 +3)*L +l];
                                const double yz = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ + l_)*I +i)*J +j)*6 +4)*L +l];
                                const double zz = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ + l_)*I +i)*J +j)*6 +5)*L +l];

                                const int xy_5d   = ((((((i_*J_ +j_)*K_ +k_)*L_ + l_)*I +i)*J +j)*5 +0)*L +l;
                                const int xz_5d   = ((((((i_*J_ +j_)*K_ +k_)*L_ + l_)*I +i)*J +j)*5 +1)*L +l;
                                const int yz_5d   = ((((((i_*J_ +j_)*K_ +k_)*L_ + l_)*I +i)*J +j)*5 +2)*L +l;
                                const int xxyy_5d = ((((((i_*J_ +j_)*K_ +k_)*L_ + l_)*I +i)*J +j)*5 +3)*L +l;
                                const int rr_5d   = ((((((i_*J_ +j_)*K_ +k_)*L_ + l_)*I +i)*J +j)*5 +4)*L +l;
                                pOutput[xy_5d] = xy;
                                pOutput[xz_5d] = xz;
                                pOutput[yz_5d] = yz;
                                pOutput[xxyy_5d] = 0.5*(xx-yy);
                                pOutput[rr_5d] = INV_SQRT3 * (zz - 0.5 * (xx + yy));
                            }
                        }
                    }
                }
            }
        }
    }
}


void DfEriEngine::transform6Dto5D_l(const int I_, const int J_, const int K_, const int L_,
                              const int I, const int J, const int K,
                              const double* pInput, double* pOutput)
{
#ifdef DEBUG_TRANSFORM_6D_TO_5D
    std::cerr << TlUtils::format("DfEriEngine::transform6Dto5D_l() in: I^=%d, J^=%d, K^=%d, L^=%d, I=%d, J=%d, K=%d",
                                 I_, J_, K_, L_, I, J, K)
              << std::endl;
#endif // DEBUG_TRANSFORM_6D_TO_5D

    for (int i_ = 0; i_ < I_; ++i_) {
        for (int j_ = 0; j_ < J_; ++j_) {
            for (int k_ = 0; k_ < K_; ++k_) {
                for (int l_ = 0; l_ < L_; ++l_) {

                    for (int i = 0; i < I; ++i) {
                        for (int j = 0; j < J; ++j) {
                            for (int k = 0; k < K; ++k) {
                                
                                const double xx = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ + l_)*I +i)*J +j)*K +k)*6 +0];
                                const double xy = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ + l_)*I +i)*J +j)*K +k)*6 +1];
                                const double xz = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ + l_)*I +i)*J +j)*K +k)*6 +2];
                                const double yy = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ + l_)*I +i)*J +j)*K +k)*6 +3];
                                const double yz = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ + l_)*I +i)*J +j)*K +k)*6 +4];
                                const double zz = pInput[((((((i_*J_ +j_)*K_ +k_)*L_ + l_)*I +i)*J +j)*K +k)*6 +5];

                                const int xy_5d   = ((((((i_*J_ +j_)*K_ +k_)*L_ + l_)*I +i)*J +j)*K +k)*5 +0;
                                const int xz_5d   = ((((((i_*J_ +j_)*K_ +k_)*L_ + l_)*I +i)*J +j)*K +k)*5 +1;
                                const int yz_5d   = ((((((i_*J_ +j_)*K_ +k_)*L_ + l_)*I +i)*J +j)*K +k)*5 +2;
                                const int xxyy_5d = ((((((i_*J_ +j_)*K_ +k_)*L_ + l_)*I +i)*J +j)*K +k)*5 +3;
                                const int rr_5d   = ((((((i_*J_ +j_)*K_ +k_)*L_ + l_)*I +i)*J +j)*K +k)*5 +4;
                                pOutput[xy_5d] = xy;
                                pOutput[xz_5d] = xz;
                                pOutput[yz_5d] = yz;
                                pOutput[xxyy_5d] = 0.5*(xx-yy);
                                pOutput[rr_5d] = INV_SQRT3 * (zz - 0.5 * (xx + yy));
                            }
                        }
                    }
                }
            }
        }
    }
}


void DfEriEngine::compD(const DfEriEngine::Query& qAB,
                        const DfEriEngine::Query& qCD)
{
    const int a = qAB.a;
    const int b = qAB.b;
    const int c = qCD.a;
    const int d = qCD.b;
    
    const int I = 2 * a + 1;
    const int J = 2 * b + 1;
    const int K = 2 * c + 1;
    const int L = 2 * d + 1;

    const int maxSize = 3 * I * J * K * L;
    for (int i = 0; i < maxSize; ++i) {
        this->WORK_D[i] = - (this->WORK_A[i] + this->WORK_B[i] + this->WORK_C[i]);
    }
}


void DfEriEngine::calc(const DfEriEngine::Query& qAB,
                       const DfEriEngine::Query& qCD)
{
    // this->log_.debug(TlUtils::format("DfEriEngine::calc(): a~=%d, b~=%d, a=%d, b=%d,",
    //                                  qAB.a_bar, qAB.b_bar, qAB.a, qAB.b));
    // this->log_.debug(TlUtils::format("                     c~=%d, d~=%d, c=%d, d=%d,",
    //                                  qCD.a_bar, qCD.b_bar, qCD.a, qCD.b));
    
    this->calcE4CQ();
    this->calc0m();
    this->calcR0();

    // choice
    const ContractScalesVector bra_contractScales_vtr = this->choice(qAB);
    const ContractScalesVector ket_contractScales_vtr = this->choice(qCD);
    
#ifdef DEBUG_CHOICE
    // for debug
    {
        int count = 0;
        for (ContractScalesType::const_iterator p = bra_contractScales_.begin();
             p != bra_contractScales_.end(); ++p) {
            std::cerr << TlUtils::format("choice[%d]: (a', b', p')=(%d, %d, %d)",
                                         count, p->a_prime, p->b_prime, p->p_prime)
                      << std::endl;
            ++count;
        }
        count = 0;
        for (ContractScalesType::const_iterator p = ket_contractScales_.begin();
             p != ket_contractScales_.end(); ++p) {
            std::cerr << TlUtils::format("choice[%d]: (c', d', q')=(%d, %d, %d)",
                                         count, p->a_prime, p->b_prime, p->p_prime)
                      << std::endl;
            ++count;
        }
    }
#endif // CHOICE
    
    // contract
    this->contract(qAB, qCD,
                   bra_contractScales_vtr,
                   ket_contractScales_vtr);

    // calc PQ
    this->calcPQ(qAB, qCD,
                 bra_contractScales_vtr,
                 ket_contractScales_vtr);

    // calc ERI for bra-
    this->calcERI(qAB.a_bar, qAB.b_bar,
                  qAB.a,     qAB.b,     0,
                  0, 0, 0,
                  &(this->ERI_bra_));

#ifdef DEBUG_BRA_ERI
    // for debug
    const ERI_State bra_state(a_bar, b_bar, a, b, 0, 0, 0, 0);
    const TlAngularMomentumVectorSet amvs_a_bar(a_bar);
    const TlAngularMomentumVectorSet amvs_b_bar(b_bar);
    const TlAngularMomentumVectorSet amvs_a(a);
    const TlAngularMomentumVectorSet amvs_b(b);
    const TlAngularMomentumVector amv_p(0, 0, 0);

    ContractScalesType::const_iterator itEnd = this->bra_contractScales_.end();
    for (ContractScalesType::const_iterator it = this->bra_contractScales_.begin(); it != itEnd; ++it) {
        const int a_prime = it->pack.a_prime;
        const int b_prime = it->pack.b_prime;
        const int p_prime = it->pack.p_prime;
        
        for (int index_a_bar = 0; index_a_bar < amvs_a_bar.size(); ++index_a_bar) {
            const TlAngularMomentumVector amv_a_bar = amvs_a_bar.get(index_a_bar);

            for (int index_b_bar = 0; index_b_bar < amvs_b_bar.size(); ++index_b_bar) {
                const TlAngularMomentumVector amv_b_bar = amvs_b_bar.get(index_b_bar);

                for (int index_a = 0; index_a < amvs_a.size(); ++index_a) {
                    const TlAngularMomentumVector amv_a = amvs_a.get(index_a);

                    for (int index_b = 0; index_b < amvs_b.size(); ++index_b) {
                        const TlAngularMomentumVector amv_b = amvs_b.get(index_b);

                        const int bra_index = this->index(amv_a_bar, amv_b_bar,
                                                          amv_a, amv_b, amv_p);
                        const std::string str = TlUtils::format("BRA_ERI: (%s %s, %s %s %s, %d %d %d|",
                                                                amv_a_bar.debugOut().c_str(),
                                                                amv_b_bar.debugOut().c_str(),
                                                                amv_a.debugOut().c_str(),
                                                                amv_b.debugOut().c_str(),
                                                                amv_p.debugOut().c_str(),
                                                                a_prime, b_prime, p_prime);
                        const std::vector<double> values = eri[bra_state][bra_index];
                        for (int i = 0; i < values.size(); ++i) {
                            std::cerr << TlUtils::format("%s [%3d] = % f",
                                                         str.c_str(), i, values[i])
                                      << std::endl;
                        }
                    }
                }
            }
        }
    }
#endif // DEBUG_BRA_ERI
    
    // swap
    this->transpose(qAB, qCD, ket_contractScales_vtr);
    std::swap(this->AB_, this->CD_); 

    // calc ERI for ket-
    this->calcERI(qCD.a_bar, qCD.b_bar,
                  qCD.a,     qCD.b,     0,
                  0, 0, 0,
                  &(this->ERI_ket_));
    std::swap(this->AB_, this->CD_); // AB, CDを元に戻す(batchブロック再計算のため)
}


void DfEriEngine::calcGrad(const DfEriEngine::Query& qAB,
                           const DfEriEngine::Query& qCD)
{
    const Query qAB10(1, 0, qAB.a, qAB.b);
    const Query qAB01(0, 1, qAB.a, qAB.b);
    const Query qAB00(0, 0, qAB.a, qAB.b);
    const Query qCD10(1, 0, qCD.a, qCD.b);
    const Query qCD00(0, 0, qCD.a, qCD.b);

    this->calcE4CQ();
    this->calc0m();
    this->calcR0();


    // AB
    const ContractScalesVector bra_contractScales_vtr = this->choice(qAB10, qAB01);
    const ContractScalesVector ket_contractScales_vtr = this->choice(qCD00);
    
    // contract
    this->contract(qAB10, qCD00,
                   bra_contractScales_vtr,
                   ket_contractScales_vtr);
    
    // A ----------------------------------------------------
    {
        // calc PQ
        this->calcPQ(qAB10, qCD00,
                     bra_contractScales_vtr,
                     ket_contractScales_vtr);
        
        // calc ERI for bra-
        this->calcERI(1, 0,
                      qAB.a,     qAB.b,     0,
                      0, 0, 0,
                      &(this->ERI_bra_));
        
        // swap
        this->transpose(qAB10, qCD00,
                        ket_contractScales_vtr);
        std::swap(this->AB_, this->CD_); 
        
        // calc ERI for ket-
        this->calcERI(0, 0,
                      qCD.a,     qCD.b,     0,
                      0, 0, 0,
                      &(this->ERI_ket_));
        std::swap(this->AB_, this->CD_); // AB, CDを元に戻す(batchブロック再計算のため)
    }
    this->copyResultsToOutputBuffer(qAB10, qCD00, this->WORK_A);
    
    // B ----------------------------------------------------
    {
        // // calc PQ
        this->calcPQ(qAB01, qCD00,
                     bra_contractScales_vtr,
                     ket_contractScales_vtr);
        
        // calc ERI for bra-
        this->calcERI(0, 1,
                      qAB.a,     qAB.b,     0,
                      0, 0, 0,
                      &(this->ERI_bra_));
        
        // swap
        this->transpose(qAB01, qCD00,
                        ket_contractScales_vtr);
        std::swap(this->AB_, this->CD_); 
        
        // calc ERI for ket-
        this->calcERI(0, 0,
                      qCD.a,     qCD.b,     0,
                      0, 0, 0,
                      &(this->ERI_ket_));
        std::swap(this->AB_, this->CD_); // AB, CDを元に戻す(batchブロック再計算のため)
    }
    this->copyResultsToOutputBuffer(qAB01, qCD00, this->WORK_B);

    // C -----------------------------------------------------
    this->calcGrad_sub(qAB00, qCD10);
    this->copyResultsToOutputBuffer(qAB00, qCD10, this->WORK_C);
}


void DfEriEngine::calcGrad_sub(const DfEriEngine::Query& qAB,
                               const DfEriEngine::Query& qCD)
{
    // choice
    const ContractScalesVector bra_contractScales_vtr = this->choice(qAB);
    const ContractScalesVector ket_contractScales_vtr = this->choice(qCD);
    
    // contract
    this->contract(qAB, qCD,
                   bra_contractScales_vtr,
                   ket_contractScales_vtr);

    // calc PQ
    this->calcPQ(qAB, qCD,
                 bra_contractScales_vtr,
                 ket_contractScales_vtr);

    // calc ERI for bra-
    this->calcERI(qAB.a_bar, qAB.b_bar,
                  qAB.a,     qAB.b,     0,
                  0, 0, 0,
                  &(this->ERI_bra_));

    // swap
    this->transpose(qAB, qCD,
                    ket_contractScales_vtr);
    std::swap(this->AB_, this->CD_); 

    // calc ERI for ket-
    this->calcERI(qCD.a_bar, qCD.b_bar,
                  qCD.a,     qCD.b,     0,
                  0, 0, 0,
                  &(this->ERI_ket_));
    std::swap(this->AB_, this->CD_); // AB, CDを元に戻す(batchブロック再計算のため)

}


// elementary 4-center quantities(E4CQ_)を作成する
// KPKQはKPが内側。
void DfEriEngine::calcE4CQ()
{
    int KP = this->bra_.size();
    int KQ = this->ket_.size();
    assert((KP * KQ) < ERI_KPKQ_MAX);

    // cutoff
    // this->bra_, this->ket_オブジェクトは大きい順に並んでいるので
    // (see DfEriEngine::getCGTO_pair())
    // 1番目(最大)に対して比較すれば良い。
    // {
    //     const double maxU_P = this->bra_[0].U_P;
    //     const double thresholdForQ = this->primitiveLevelThreshold_ / std::fabs(maxU_P);
    //     for (int i = 0; i < KQ; ++i) {
    //         const double U_Q = this->ket_[i].U_P;
    //         if (std::fabs(U_Q) <= thresholdForQ) {
    //             KQ = i + 1;
    //             this->ket_.resize(KQ);
    //             break;
    //         }
    //     }
    //     const double maxU_Q = this->ket_[0].U_P;
    //     const double thresholdForP = this->primitiveLevelThreshold_ / std::fabs(maxU_Q);
    //     for (int i = 0; i < KP; ++i) {
    //         const double U_P = this->bra_[i].U_P;
    //         if (std::fabs(U_P) <= thresholdForP) {
    //             KP = i + 1;
    //             this->bra_.resize(KP);
    //             break;
    //         }
    //     }
    // }
    

    // calc
    int index = 0;
    for (int KQ_index = 0; KQ_index < KQ; ++KQ_index) {
        const TlPosition Q = this->ket_[KQ_index].P();
        const double sigma_Q = this->ket_[KQ_index].sigma_P();
        const double U_Q = this->ket_[KQ_index].U_P();

        for (int KP_index = 0; KP_index < KP; ++KP_index) {
            const TlPosition P = this->bra_[KP_index].P();
            const double sigma_P = this->bra_[KP_index].sigma_P();
            const double U_P = this->bra_[KP_index].U_P();

            assert(index == KP * KQ_index + KP_index);
            this->pE4_[index] = E4(P, Q, sigma_P, sigma_Q, U_P, U_Q);
            ++index;
        }
    }
}


// [0]^(m) を求める
// 計算結果の格納先は _0M_
void DfEriEngine::calc0m()
{
    const int sumOfAngularMomentums =  this->sumOfAngularMomentums_;
    assert(sumOfAngularMomentums < ERI_L_MAX);
    
    // Gm(T)を求める
    const TlFmt& fmt = TlFmt::getInstance();
    const int KPKQ = this->bra_.size() * this->ket_.size();
    for (int kpkq = 0; kpkq < KPKQ; ++kpkq) {
        const double T = 0.5 * this->pE4_[kpkq]._2T;
        fmt.getGmT(sumOfAngularMomentums, T, this->pGmT_);

        const double U = this->pE4_[kpkq].U;
        const double _2theta2 = this->pE4_[kpkq]._2theta2;
        
        double pow_2theta2 = std::sqrt(_2theta2); // = pow(_2theta2, m+1/2)
        for (int m = 0; m <= sumOfAngularMomentums; ++m) {
            this->p0M_[m][kpkq] = U * pow_2theta2 * this->pGmT_[m];
            pow_2theta2 *= _2theta2;
        }
    }
}


void DfEriEngine::calcR0()
{
    const int M = this->sumOfAngularMomentums_;
    const int KPKQ = this->bra_.size() * this->ket_.size();
    
    // 初期化 ------------------------------------------------------------------
    this->calcdRM_.reset();
    
    // p0M から値をコピー ------------------------------------------------------
    const TlAngularMomentumVector r(0, 0, 0);
    for (int m = 0; m <= M; ++m) {
        const unsigned int index = this->indexRM(r, m);
        std::copy(this->p0M_[m], this->p0M_[m] + KPKQ, this->pRM_[index]);
        this->calcdRM_[index] = true;
    }

#ifdef DEBUG_CALC_R0
    // for debug
    {
        const TlAngularMomentumVector zero(0, 0, 0);
        for (int m = 0; m <= M; ++m) {
            const unsigned int index = this->indexRM(zero, m);
            for (int i = 0; i < KPKQ; ++i) {
                std::cerr << TlUtils::format("calcR0(in): [0](%d) [%d] = %f",
                                             m, i, this->pRM_[index][i])
                          << std::endl;
            }
        }
    }
#endif // DEBUG_CALC_R0
    
    // [r]^(0)を計算 -----------------------------------------------------------
    for (int am = 0; am <= M; ++am) { // am means Angular Momentum
        const TlAngularMomentumVectorSet amvs(am);
        const int numOfAMVS = amvs.size();
        for (int i = 0; i < numOfAMVS; ++i) {
            const TlAngularMomentumVector amv = amvs.get(i);
            this->calcRM(amv, 0); // [r]^(0)を計算
        }
    }
    
#ifdef DEBUG_CALC_R0
    // for debug 
    for (int am = 0; am <= M; ++am) {
        const TlAngularMomentumVectorSet amvs(am);
        int size = amvs.size();
        for (int i = 0; i < size; ++i) {
            const TlAngularMomentumVector amv = amvs.get(i);
            const int index = this->indexRM(amv, 0);
            for (int kpkq = 0; kpkq < KPKQ; ++kpkq) {
                std::cerr << TlUtils::format("calcR(out): %s(0) [%d] = %f",
                                             amv.debugOut().c_str(), kpkq, this->pRM_[index][kpkq])
                          << std::endl;
            }
        }
    }
#endif // DEBUG_CALC_R0
}


unsigned int DfEriEngine::indexRM(const TlAngularMomentumVector& amv, const int m) const
{
    const int L = amv.angularMomentum();
    assert((L + m) < ERI_L_MAX);
    const unsigned int amvBegin = L*(L+1)*(L+2)/6; // [0](m)から[L-1](m)までが何個あったか
    const unsigned int index = ERI_NUM_OF_R_KINDS * m + amvBegin + amv.index();
    
    return index;
}


void DfEriEngine::calcRM(const TlAngularMomentumVector& r, const int m) {
    const unsigned int index = this->indexRM(r, m);

    if (this->calcdRM_[index] != true) {
        if (r.angularMomentum() > 0) {
            const int KPKQ = this->bra_.size() * this->ket_.size();
            const int i = this->initiativeRM(r);
            // 1st term
            {
                const TlAngularMomentumVector r1 = r - this->E1_[i];
                this->calcRM(r1, m + 1);
                const unsigned int index1 = this->indexRM(r1, m + 1);
                
                for (int k = 0; k < KPKQ; ++k) {
                    const double coef = this->pE4_[k].R[i];
                    this->pRM_[index][k] = coef * this->pRM_[index1][k];
                }
            }
            
            // 2nd term
            const TlAngularMomentumVector r2 = r - this->E2_[i];
            if (r2.isExist() == true) {
                this->calcRM(r2, m + 1);
                const unsigned int index2 = this->indexRM(r2, m + 1);
                
                const double coef = r.get(i) - 1.0;
                for (int k = 0; k < KPKQ; ++k) {
                    this->pRM_[index][k] -= coef * this->pRM_[index2][k];
                }
            }
        }
    
        this->calcdRM_[index] = true;
    }
}


int DfEriEngine::initiativeRM(const TlAngularMomentumVector& amv) const {
    //const int x = amv.get(0);
    const int y = amv.get(1);
    const int z = amv.get(2);
    
    if (z > 0) {
        return 2;
    } else if (y > 0) {
        return 1;
    } else {
        return 0;
    }
}

DfEriEngine::ContractScalesVector
DfEriEngine::choice(const DfEriEngine::Query& qAB)
{
#ifdef USE_CACHED_ROUTE
    const int index = qAB.index();
    if (this->choice_tbl_.find(index) == this->choice_tbl_.end()) {
        ContractScalesSet contractList;
        this->choice(qAB.a_bar, qAB.b_bar,
                     qAB.a, qAB.b, 0,
                     0, 0, 0,
                     &contractList);
        this->choice_tbl_[index] = this->transContractScales_SetToVector(contractList);
    }
    return this->choice_tbl_[index];
#else
    {
        ContractScalesSet contractList;
        this->choice(qAB.a_bar, qAB.b_bar,
                     qAB.a, qAB.b, 0,
                     0, 0, 0,
                     &contractList);
        return this->transContractScales_SetToVector(contractList);
    }
#endif // USE_CACHED_ROUTE
}


DfEriEngine::ContractScalesVector
DfEriEngine::choice(const DfEriEngine::Query& qAB1,
                    const DfEriEngine::Query& qAB2)
{
    ContractScalesSet contractList;
    this->choice(qAB1.a_bar, qAB1.b_bar,
                 qAB1.a, qAB1.b, 0,
                 0, 0, 0,
                 &contractList);
    this->choice(qAB2.a_bar, qAB2.b_bar,
                 qAB2.a, qAB2.b, 0,
                 0, 0, 0,
                 &contractList);
    return this->transContractScales_SetToVector(contractList);
}


// pContractList に計算に必要な状態が登録される
void DfEriEngine::choice(const int a_bar, const int b_bar,
                         const int a, const int b, const int p,
                         const int a_prime, const int b_prime, const int p_prime,
                         ContractScalesSet* pContractList)
{
    assert(pContractList != NULL);

    if ((a_bar < 0) || (b_bar < 0) ||
        (a < 0) || (b < 0) || (p < 0) ||
        (a_prime < 0) || (b_prime < 0) || (p_prime < 0)) {
        // do nothing because this state does not exist!
        return;
    }

#ifdef CHECK_MAX_COUNT
    this->max_a_bar_ = std::max(this->max_a_bar_, a_bar);
    this->max_b_bar_ = std::max(this->max_b_bar_, b_bar);
    this->max_a_ = std::max(this->max_a_, a);
    this->max_b_ = std::max(this->max_b_, b);
    this->max_p_ = std::max(this->max_p_, p);
    this->max_a_prime_ = std::max(this->max_a_prime_, a_prime);
    this->max_b_prime_ = std::max(this->max_b_prime_, b_prime);
    this->max_p_prime_ = std::max(this->max_p_prime_, p_prime);
#endif // CHECK_MAX_COUNT

    // recursive
    if ((a_bar == 0) && (b_bar == 0) && (a == 0) && (b == 0)) {
        // register for contract
        pContractList->insert(ContractScale(a_prime, b_prime, p_prime));

        // stop
        return;
    } else if (b_bar > 0) {
        // use eq.47
        this->choice(a_bar,   b_bar-1, a,   b,   p+1, a_prime,   b_prime,   p_prime,   pContractList);
        this->choice(a_bar+1, b_bar-1, a,   b,   p,   a_prime,   b_prime,   p_prime,   pContractList);
        this->choice(a_bar,   b_bar-1, a-1, b,   p,   a_prime,   b_prime,   p_prime,   pContractList);
        this->choice(a_bar,   b_bar-1, a,   b-1, p,   a_prime,   b_prime,   p_prime,   pContractList);
    } else if (b > 0) {
        // use eq.44
        this->choice(a_bar,   b_bar,   a+1, b-1, p,   a_prime,   b_prime,   p_prime,   pContractList);
        this->choice(a_bar,   b_bar,   a,   b-1, p,   a_prime,   b_prime,   p_prime,   pContractList);
        this->choice(a_bar-1, b_bar,   a,   b-1, p,   a_prime,   b_prime,   p_prime,   pContractList);
        this->choice(a_bar,   b_bar-1, a,   b-1, p,   a_prime,   b_prime,   p_prime,   pContractList);
    } else if (a_bar > 0) {
        // use eq.43
        this->choice(a_bar-1, b_bar,   a+1, b,   p,   a_prime+1, b_prime,   p_prime,   pContractList);
        this->choice(a_bar-1, b_bar,   a-1, b,   p,   a_prime,   b_prime,   p_prime,   pContractList);
        this->choice(a_bar-1, b_bar,   a,   b,   p-1, a_prime+1, b_prime,   p_prime,   pContractList);
    } else if (a > 0) {
        // use eq.45 (p -> a)
        this->choice(a_bar,   b_bar,   a-1, b,   p-1, a_prime,   b_prime,   p_prime  , pContractList);
        this->choice(a_bar,   b_bar,   a-1, b,   p+1, a_prime,   b_prime,   p_prime+1, pContractList);
        this->choice(a_bar,   b_bar,   a-1, b,   p,   a_prime,   b_prime+1, p_prime+1, pContractList);
        this->choice(a_bar-1, b_bar,   a-1, b,   p,   a_prime,   b_prime+1, p_prime+1, pContractList);
        this->choice(a_bar,   b_bar-1, a-1, b,   p,   a_prime,   b_prime+1, p_prime+1, pContractList);
    } else {
        // something wrong
        std::cerr << TlUtils::format("program error. @ DfEriEngine::choice() [%d, %d, %d, %d, %d, %d, %d, %d]",
                                     a_bar, b_bar, a, b, p, a_prime, b_prime, p_prime)
                  << std::endl;
    }
}


DfEriEngine::ContractScalesVector
DfEriEngine::transContractScales_SetToVector(const ContractScalesSet& contractScales)
{
    const std::size_t size = contractScales.size();
    ContractScalesVector answer(size);

    std::size_t index = 0;
    ContractScalesSet::const_iterator pEnd = contractScales.end();
    for (ContractScalesSet::const_iterator p = contractScales.begin();  p != pEnd; ++p) {
        answer[index] = *p;
        ++index;
    }

    return answer;
}


void DfEriEngine::contract(const DfEriEngine::Query& qAB,
                           const DfEriEngine::Query& qCD,
                           const ContractScalesVector& bra_contractScales,
                           const ContractScalesVector& ket_contractScales)
{
    const int R = this->sumOfAngularMomentums_;
    //this->nR_dash_.clear();
    
    // contract bra-
    std::size_t nR_dash_index = 0;
    const int max_bra_cs_index = bra_contractScales.size();
    for (int bra_cs_index = 0; bra_cs_index < max_bra_cs_index; ++bra_cs_index) {
        const int a_prime = bra_contractScales[bra_cs_index].a_prime;
        const int b_prime = bra_contractScales[bra_cs_index].b_prime;
        const int p_prime = bra_contractScales[bra_cs_index].p_prime;
        
        for (int r = 0; r <= R; ++r) {
//             std::cerr << TlUtils::format("contract() for bra- a'=%d b'=%d p'=%d r=%d",
//                                          a_prime, b_prime, p_prime, r) << std::endl;
            const TlAngularMomentumVectorSet amvs(r);
            const int numOfAmv = amvs.size();
            for (int amv_index = 0; amv_index < numOfAmv; ++amv_index) {
                const TlAngularMomentumVector amv = amvs.get(amv_index);
                
                assert(nR_dash_index < ERI_NR_DASH_SIZE);
                this->contract_bra(qAB, amv, a_prime, b_prime, p_prime, nR_dash_index);
                ++nR_dash_index;
            }
        }
    }
    
    // for profile
#ifdef CHECK_MAX_COUNT
    this->maxSizeOf_nR_dash_ = std::max(this->maxSizeOf_nR_dash_, nR_dash_index);
#endif // SHOW_MAX_COUNT
    
    // for debug
#ifdef CHECK_R_DASH
    for (n_R_dash_Type::const_iterator p = this->n_R_dash_.begin();
         p != this->n_R_dash_.end(); ++p) {
        ContractState cs = p->first;
        const std::vector<double> values = p->second;
        std::string v_string = "[";
        for (std::vector<double>::const_iterator q = values.begin(); q != values.end(); ++q) {
            v_string += TlUtils::format("%e ", *q);
        }
        v_string += "]";
        std::cerr << TlUtils::format("(R') %s = %s",
                                     cs.debugOut().c_str(),
                                     v_string.c_str())
                  << std::endl;
    }
#endif //DEBUG_R_DASH

    // contract ket-
    const int cd = qCD.a + qCD.b;
    const int KQ = this->ket_.size();
    //assert(KQ == (int)KQ_values.size());
    std::vector<double> coef_numerators(KQ);
    const int max_ket_cs_index = ket_contractScales.size();
    for (int ket_cs_index = 0; ket_cs_index < max_ket_cs_index; ++ket_cs_index) {
        const int c_prime = ket_contractScales[ket_cs_index].a_prime;
        const int d_prime = ket_contractScales[ket_cs_index].b_prime;
        const int q_prime = ket_contractScales[ket_cs_index].p_prime;

        const int zeta_exp = q_prime - cd;
        this->get_contract_ket_coef_numerators(c_prime, d_prime, zeta_exp,
                                               &coef_numerators);

        for (std::size_t i = 0; i < nR_dash_index; ++i) {
            ContractState cs = this->nR_dash_[i].cs;
            cs.setCDQ(c_prime, d_prime, q_prime);
            const std::vector<double>& values = this->nR_dash_[i].values;

            const std::size_t cs_index = cs.index();
            assert(cs_index < (ERI_P_PRIME_MAX * ERI_P_PRIME_MAX * ERI_P_PRIME_MAX *
                               ERI_P_PRIME_MAX * ERI_P_PRIME_MAX * ERI_P_PRIME_MAX * 
                               ERI_NUM_OF_R_KINDS));

            double value = 0.0;
            for (int KQ_index = 0; KQ_index < KQ; ++KQ_index) {
                value += coef_numerators[KQ_index] * values[KQ_index];
            }

            this->p_abpRcdq_[cs_index] = value;
        }
    }

//     // for debug
// #ifdef DEBUG_R_2DASH
//     for (n_R_dash2_Type::const_iterator p = this->n_R_dash2_.begin();
//          p != this->n_R_dash2_.end(); ++p) {
//         ContractState cs = p->first;
//         double value = p->second;
//         std::cerr << TlUtils::format("(R'') %s = % f",
//                                      cs.debugOut().c_str(),
//                                      value)
//                   << std::endl;
//     }
// #endif // DEBUG_R_2DASH
}


// eq.36
void DfEriEngine::contract_bra(const DfEriEngine::Query& qAB,
                               const TlAngularMomentumVector& r,
                               const int a_prime, const int b_prime, const int p_prime,
                               const std::size_t nR_dash_index)
{
    const int a = qAB.a;
    const int b = qAB.b;
    const int zeta_exp = p_prime - a - b;

    const int KP = this->bra_.size();
    const int KQ = this->ket_.size();
    const unsigned int indexR0 = this->indexRM(r, 0);

    //this->pContractBraCoef_ はメモリ確保の回数を減らすためにメンバ変数へ
    assert(KP < ERI_KPKQ_MAX);
    assert(KQ < ERI_KPKQ_MAX);
    
    ContractState& cs = this->nR_dash_[nR_dash_index].cs;
    // cs.r_index = r.index();
    // cs.r = r.angularMomentum();
    cs.setR(r);
    // cs.a_prime = a_prime;
    // cs.b_prime = b_prime;
    // cs.p_prime = p_prime;
    cs.setABP(a_prime, b_prime, p_prime);
    std::vector<double>& values = this->nR_dash_[nR_dash_index].values;
    values.resize(KQ);

    for (int KP_index = 0; KP_index < KP; ++KP_index) {
        const double alpha2 = this->bra_[KP_index].alpha2();
        const double beta2  = this->bra_[KP_index].beta2();
        const double zeta2  = this->bra_[KP_index].zeta2();
        const double _2a   = TlMath::pow(alpha2, a_prime);
        const double _2b   = TlMath::pow(beta2,  b_prime);

        const double _2z   = TlMath::pow(zeta2,  - zeta_exp);
        this->pContractBraCoef_[KP_index] = _2a * _2b * _2z;
        // const double _2z   = TlMath::pow(zeta2,  zeta_exp);
        // this->pContractBraCoef_[KP_index] = _2a * _2b / _2z;
    }

    for (int KQ_index = 0; KQ_index < KQ; ++KQ_index) {

        double value = 0.0;
        for (int KP_index = 0; KP_index < KP; ++KP_index) {
            const int KPKQ_index = KP * KQ_index + KP_index;
            const double r0 = this->pRM_[indexR0][KPKQ_index];
            value += this->pContractBraCoef_[KP_index] * r0;
        }
        
#ifdef DEBUG_CONTRACT_BRA
        std::cerr << TlUtils::format("contract_bra[%2d]: alpha=% f, beta=% f, zeta=% f, exp_zeta=%d, coef=% f, r0=% f",
                                     KP_index, alpha, beta, zeta, zeta_exp, coef, r0)
                  << std::endl;
#endif // DEBUG_CONTRACT_BRA
        
        values[KQ_index] = value;
    }
    
#ifdef DEBUG_CONTRACT_BRA
    for (int i = 0; i < tmp.size(); ++i) {
        std::cerr << TlUtils::format("(R') %s [%d] = % f",
                                     cs.debugOut().c_str(), i, values[i])
                  << std::endl;
    }
#endif // DEBUG_CONTRACT_BRA
}

// eq.36
void DfEriEngine::get_contract_ket_coef_numerators(const int c_prime,
                                                   const int d_prime,
                                                   const int zeta_exp,
                                                   std::vector<double>* pCoefNumerators)
{
    const int KQ = this->ket_.size();
    for (int KQ_index = 0; KQ_index < KQ; ++KQ_index) {
        const double alpha2 = this->ket_[KQ_index].alpha2();
        const double beta2 = this->ket_[KQ_index].beta2();
        const double zeta2 = this->ket_[KQ_index].zeta2();
        const double _2a = TlMath::pow(alpha2, c_prime);
        const double _2b = TlMath::pow(beta2,  d_prime);
        const double _2z = TlMath::pow(zeta2,  zeta_exp);
        (*pCoefNumerators)[KQ_index] = _2a * _2b / _2z;
    }
}

// void DfEriEngine::contract_ket(const ContractState& cs,
//                                const std::vector<double>& coef_numerators,
//                                const std::vector<double>& KQ_values)
// {
//     double value = 0.0;
//     const int KQ = this->ket_.size();
//     assert(KQ == (int)KQ_values.size());
//     for (int KQ_index = 0; KQ_index < KQ; ++KQ_index) {
//         value += coef_numerators[KQ_index] * KQ_values[KQ_index];
//     }

//     const std::size_t cs_index = cs.index();
//     assert(cs_index < (ERI_P_PRIME_MAX * ERI_P_PRIME_MAX * ERI_P_PRIME_MAX *
//                        ERI_P_PRIME_MAX * ERI_P_PRIME_MAX * ERI_P_PRIME_MAX * ERI_NUM_OF_R_KINDS));
//     this->p_abpRcdq_[cs_index] = value;
// }

// eq.36
// void DfEriEngine::contract_ket(const DfEriEngine::Query& qCD,
//                                const ContractState& cs, const std::vector<double>& KQ_values)
// {
//     const int c = qCD.a;
//     const int d = qCD.b;

//     // const int c_prime = cs.c_prime;
//     // const int d_prime = cs.d_prime;
//     // const int q_prime = cs.q_prime;
//     const int c_prime = cs.getCprime();
//     const int d_prime = cs.getDprime();
//     const int q_prime = cs.getQprime();
//     const int zeta_exp = q_prime - c - d;
    
//     double value = 0.0;
//     const int KQ = this->ket_.size();
//     assert(KQ == (int)KQ_values.size());
//     for (int KQ_index = 0; KQ_index < KQ; ++KQ_index) {
//         const double alpha2 = this->ket_[KQ_index].alpha2();
//         const double beta2 = this->ket_[KQ_index].beta2();
//         const double zeta2 = this->ket_[KQ_index].zeta2();
//         const double _2a = TlMath::pow(alpha2, c_prime);
//         const double _2b = TlMath::pow(beta2,  d_prime);
//         const double _2z = TlMath::pow(zeta2,  zeta_exp);
//         const double coef = _2a * _2b / _2z;

//         value += coef * KQ_values[KQ_index];
//     }

//     const std::size_t cs_index = cs.index();
//     assert(cs_index < (ERI_P_PRIME_MAX * ERI_P_PRIME_MAX * ERI_P_PRIME_MAX *
//                        ERI_P_PRIME_MAX * ERI_P_PRIME_MAX * ERI_P_PRIME_MAX * ERI_NUM_OF_R_KINDS));
//     this->p_abpRcdq_[cs_index] = value;
// }


void DfEriEngine::calcPQ(const DfEriEngine::Query& qAB,
                         const DfEriEngine::Query& qCD,
                         const ContractScalesVector& bra_contractScales,
                         const ContractScalesVector& ket_contractScales)
{
    this->isCalcdERI_.clear();
    
    const int angularMomentumP =
          qAB.a + qAB.b 
        + qAB.a_bar + qAB.b_bar;
    const int angularMomentumQ =
          qCD.a + qCD.b
        + qCD.a_bar + qCD.b_bar;

    const int a = 0;
    const int b = 0;
    const int a_bar = 0;
    const int b_bar = 0;
    const TlAngularMomentumVector amv_a_bar(0, 0, 0);
    const TlAngularMomentumVector amv_b_bar(0, 0, 0);
    const TlAngularMomentumVector amv_a(0, 0, 0);
    const TlAngularMomentumVector amv_b(0, 0, 0);
    const int max_ket_index = ket_contractScales.size()
        * ((angularMomentumQ +1) * (angularMomentumQ +2) * (angularMomentumQ +3) / 6);
    this->ERI_batch_ = max_ket_index;

#ifdef CHECK_MAX_COUNT
    this->maxERI_batch_ = std::max(this->maxERI_batch_, this->ERI_batch_);
#endif // CHECK_MAX_COUNT
    
#ifdef DEBUG_HGP
    std::cerr << "angularMomentumP=" << angularMomentumP << ", angularMomentumQ=" << angularMomentumQ << std::endl;
    std::cerr << "calcPQ: max_ket_index=" << max_ket_index << std::endl;
#endif // DEBUG_HGP

    // ContractState cs;
    const std::vector<int> cs_indeces = this->get_csindex_for_calcPQ(angularMomentumP, qAB,
                                                                     angularMomentumQ, qCD);
    int cs_index_order = 0;

    // bra ---------------------------------------------------------------------
    const int max_bra_cs_index = bra_contractScales.size();
    for (int bra_cs_index = 0; bra_cs_index < max_bra_cs_index; ++bra_cs_index) {
        const int a_prime = bra_contractScales[bra_cs_index].a_prime;
        const int b_prime = bra_contractScales[bra_cs_index].b_prime;
        const int p_prime = bra_contractScales[bra_cs_index].p_prime;
        // cs.setABP(a_prime, b_prime, p_prime);

        for (int p = 0; p <= angularMomentumP; ++p) {
            const ERI_State bra_state(a_bar, b_bar, a, b, p, a_prime, b_prime, p_prime);
            const std::size_t braStateIndex = bra_state.index();
            this->isCalcdERI_[bra_state] = true;

// #ifdef DEBUG_CALC_PQ
//             // for debug
//             std::cerr << TlUtils::format("calcPQ: set %s",
//                                          bra_state.debugOut().c_str())
//                       << std::endl;
// #endif // DEBUG_CALC_PQ

            const TlAngularMomentumVectorSet amvs_p(p);
            const int numOf_amvs_p = amvs_p.size(); // amvs_p の数だけで十分
            assert(braStateIndex < ERI_NUM_OF_ERI_STATES);
            this->ERI_bra_[braStateIndex].resize(numOf_amvs_p);

#ifdef CHECK_MAX_COUNT
            this->maxNumOfAMVs_ = std::max(this->maxNumOfAMVs_, numOf_amvs_p);
#endif // CHECK_MAX_COUNT
    
            for (int i = 0; i < numOf_amvs_p; ++i) {
                const TlAngularMomentumVector amv_p = amvs_p.get(i);
                const int bra_index = this->index(amv_a_bar, amv_b_bar,
                                                  amv_a, amv_b, amv_p);
                assert(bra_index < numOf_amvs_p);

                // ket ---------------------------------------------------------
                std::vector<double>& ket = this->ERI_bra_[braStateIndex][bra_index];
                ket.resize(ERI_MAX_BATCH);
                int ket_index = 0;
                const int max_ket_cs_index = ket_contractScales.size();
                for (int ket_cs_index = 0; ket_cs_index < max_ket_cs_index; ++ket_cs_index) {
                    // const int c_prime = ket_contractScales[ket_cs_index].a_prime;
                    // const int d_prime = ket_contractScales[ket_cs_index].b_prime;
                    // const int q_prime = ket_contractScales[ket_cs_index].p_prime;
                    // cs.setCDQ(c_prime, d_prime, q_prime);

                    double coef = -1.0;
                    for (int q = 0; q <= angularMomentumQ; ++q) {
                        const TlAngularMomentumVectorSet amvs_q(q);
                        //const double coef = ((q % 2) == 0) ? 1.0 : -1.0;
                        coef *= -1.0; // [1.0(q=0), -1.0(q=1), 1.0(q=2), ... ]

                        const int numOf_amvs_q = amvs_q.size();
                        for (int j = 0; j < numOf_amvs_q; ++j) {
                            // const TlAngularMomentumVector amv_q = amvs_q.get(j);
                            // const TlAngularMomentumVector r = amv_p + amv_q;
                            // cs.setR(r);

                            // const int cs_index = cs.index();
                            // if (cs_indeces[cs_index_order] != cs_index) {
                            //     std::cerr << TlUtils::format("err! %d != %d",
                            //                                  cs_indeces[cs_index_order], cs_index)
                            //               << std::endl;
                            // }
                            const int cs_index = cs_indeces[cs_index_order];
                            ++cs_index_order;

                            const double value = this->p_abpRcdq_[cs_index];
#ifdef DEBUG_CALC_PQ
                            // for debug
                            std::cerr << TlUtils::format("calcPQ: %s(p=%s) [ket:%d] = (R'') (%s: %d %d %d|%d %d %d)=% f*% f",
                                                         bra_state.debugOut().c_str(),
                                                         amv_p.debugOut().c_str(),
                                                         ket_index,
                                                         r.debugOut().c_str(),
                                                         a_prime, b_prime, p_prime,
                                                         c_prime, d_prime, q_prime,
                                                         value, coef)
                                      << std::endl;
#endif // DEBUG_CALC_PQ
                            assert(ket_index < max_ket_index);
                            ket[ket_index] = coef * value;
                            ++ket_index;
                        }
                    }
                }
                assert(max_ket_index == ket_index);
            }
        }
    }
}


std::vector<int> 
DfEriEngine::get_csindex_for_calcPQ(const int angularMomentumP,
                                    const DfEriEngine::Query& qAB,
                                    const int angularMomentumQ,
                                    const DfEriEngine::Query& qCD)
{
    const int index = 
        ((angularMomentumP*Query::maxIndex() + qAB.index())*ERI_P_MAX + angularMomentumQ)*Query::maxIndex() + qCD.index();
    if (this->csindeces_forCalcPQ_.find(index) == this->csindeces_forCalcPQ_.end()) {
        std::vector<int> answer;
        const ContractScalesVector bra_contractScales = this->choice(qAB);
        const ContractScalesVector ket_contractScales = this->choice(qCD);
        
        ContractState cs;
        // bra ---------------------------------------------------------------------
        const int max_bra_cs_index = bra_contractScales.size();
        for (int bra_cs_index = 0; bra_cs_index < max_bra_cs_index; ++bra_cs_index) {
            const int a_prime = bra_contractScales[bra_cs_index].a_prime;
            const int b_prime = bra_contractScales[bra_cs_index].b_prime;
            const int p_prime = bra_contractScales[bra_cs_index].p_prime;
            cs.setABP(a_prime, b_prime, p_prime);
            
            for (int p = 0; p <= angularMomentumP; ++p) {
                const TlAngularMomentumVectorSet amvs_p(p);
                const int numOf_amvs_p = amvs_p.size(); // amvs_p の数だけで十分
                for (int i = 0; i < numOf_amvs_p; ++i) {
                    const TlAngularMomentumVector amv_p = amvs_p.get(i);
                    
                    // ket ---------------------------------------------------------
                    int ket_index = 0;
                    const int max_ket_cs_index = ket_contractScales.size();
                    for (int ket_cs_index = 0; ket_cs_index < max_ket_cs_index; ++ket_cs_index) {
                        const int c_prime = ket_contractScales[ket_cs_index].a_prime;
                        const int d_prime = ket_contractScales[ket_cs_index].b_prime;
                        const int q_prime = ket_contractScales[ket_cs_index].p_prime;
                        cs.setCDQ(c_prime, d_prime, q_prime);
                        
                        for (int q = 0; q <= angularMomentumQ; ++q) {
                            const TlAngularMomentumVectorSet amvs_q(q);
                            const int numOf_amvs_q = amvs_q.size();
                            for (int j = 0; j < numOf_amvs_q; ++j) {
                                const TlAngularMomentumVector amv_q = amvs_q.get(j);
                                const TlAngularMomentumVector r = amv_p + amv_q;
                                cs.setR(r);
                                const int cs_index = cs.index();
                                answer.push_back(cs_index);
                            }
                        }
                    }
                }
            }
        }
        this->csindeces_forCalcPQ_[index] = answer;
    }

    return this->csindeces_forCalcPQ_[index];
}

void DfEriEngine::transpose(const DfEriEngine::Query& qAB,
                            const DfEriEngine::Query& qCD,
                            const ContractScalesVector& ket_contractScales)
{
    this->isCalcdERI_.clear();
    
    const int angularMomentumQ =
          qCD.a + qCD.b
        + qCD.a_bar + qCD.b_bar;
    
    const int a_bar = qAB.a_bar;
    const int b_bar = qAB.b_bar;
    const int a = qAB.a;
    const int b = qAB.b;
    const int p = 0;
    const int a_prime = 0;
    const int b_prime = 0;
    const int p_prime = 0;
    const ERI_State bra_state(a_bar, b_bar, a, b, p, a_prime, b_prime, p_prime);
    const std::size_t braStateIndex = bra_state.index();
    assert(braStateIndex < ERI_NUM_OF_ERI_STATES);
    
    const int c = 0;
    const int d = 0;
    const int c_bar = 0;
    const int d_bar = 0;

    const TlAngularMomentumVectorSet amvs_a_bar(a_bar);
    const TlAngularMomentumVectorSet amvs_b_bar(b_bar);
    const TlAngularMomentumVectorSet amvs_a(a);
    const TlAngularMomentumVectorSet amvs_b(b);
    const int max_amvs_a_bar_index = amvs_a_bar.size();
    const int max_amvs_b_bar_index = amvs_b_bar.size();
    const int max_amvs_a_index = amvs_a.size();
    const int max_amvs_b_index = amvs_b.size();
    const TlAngularMomentumVector amv_p(0, 0, 0);
    
    const int max_bra_index = max_amvs_a_bar_index * max_amvs_b_bar_index * max_amvs_a_index * max_amvs_b_index;
    this->ERI_batch_ = max_bra_index;

#ifdef CHECK_MAX_COUNT
    this->maxERI_batch_ = std::max(this->maxERI_batch_, this->ERI_batch_);
#endif // CHECK_MAX_COUNT

#ifdef DEBUG_HGP
    std::cerr << "DfEriEngine::transpose(): angularMomentumQ=" << angularMomentumQ << std::endl;
#endif // DEBUG_HGP
    
    // ket ---------------------------------------------------------------------
    int ket_index = 0;
    // ContractScalesVector::const_iterator it_ket_end = ket_contractScales.end();
    // for (ContractScalesVector::const_iterator it_ket = ket_contractScales.begin();
         // it_ket != it_ket_end; ++it_ket) {
    const int max_ket_cs_index = ket_contractScales.size();
    for (int ket_cs_index = 0; ket_cs_index < max_ket_cs_index; ++ket_cs_index) {
        const int c_prime = ket_contractScales[ket_cs_index].a_prime;
        const int d_prime = ket_contractScales[ket_cs_index].b_prime;
        const int q_prime = ket_contractScales[ket_cs_index].p_prime;

        for (int q = 0; q <= angularMomentumQ; ++q) {
            const ERI_State ket_state(c_bar, d_bar, c, d, q, c_prime, d_prime, q_prime);
            const std::size_t ketStateIndex = ket_state.index();
            this->isCalcdERI_[ket_state] = true;

#ifdef DEBUG_HGP
            std::cerr << TlUtils::format("DfEriEngine::transpose(): set %s", ket_state.debugOut().c_str())
                      << std::endl;
#endif // DEBUG_HGP
            
            const TlAngularMomentumVectorSet amvs_q(q);
            const int numOf_amvs_q = amvs_q.size();
            this->ERI_ket_[ketStateIndex].resize(numOf_amvs_q);

#ifdef CHECK_MAX_COUNT
            this->maxNumOfAMVs_ = std::max(this->maxNumOfAMVs_, numOf_amvs_q);
#endif // CHECK_MAX_COUNT
            
            for (int j = 0; j < numOf_amvs_q; ++j) {
                const TlAngularMomentumVector amv_q = amvs_q.get(j);
                const int amvs_q_index = amv_q.index();
                //this->ERI_ket_[ket_state.index()][amvs_q_index].resize(max_bra_index);
                this->ERI_ket_[ket_state.index()][amvs_q_index].resize(ERI_MAX_BATCH);
                
                // bra ---------------------------------------------------------
                int bra_index = 0;
                for (int amvs_a_bar_index = 0; amvs_a_bar_index < max_amvs_a_bar_index; ++amvs_a_bar_index) {
                    const TlAngularMomentumVector amv_a_bar = amvs_a_bar.get(amvs_a_bar_index);

                    for (int amvs_b_bar_index = 0; amvs_b_bar_index < max_amvs_b_bar_index; ++amvs_b_bar_index) {
                        const TlAngularMomentumVector amv_b_bar = amvs_b_bar.get(amvs_b_bar_index);

                        for (int amvs_a_index = 0; amvs_a_index < max_amvs_a_index; ++amvs_a_index) {
                            const TlAngularMomentumVector amv_a = amvs_a.get(amvs_a_index);

                            for (int amvs_b_index = 0; amvs_b_index < max_amvs_b_index; ++amvs_b_index) {
                                const TlAngularMomentumVector amv_b = amvs_b.get(amvs_b_index);

                                const int amvs_bra_index = this->index(amv_a_bar, amv_b_bar,
                                                                       amv_a, amv_b, amv_p);

//                                 std::cerr << TlUtils::format("transpose() copy from [%s %d %d] to [%s %d %d], ",
//                                                              bra_state.debugOut().c_str(), amvs_bra_index, ket_index,
//                                                              ket_state.debugOut().c_str(), amvs_ket_index, bra_index);

                                const double value = this->ERI_bra_[braStateIndex][amvs_bra_index][ket_index];

//                                 std::cerr << "value = " << value << std::endl;
                                this->ERI_ket_[ketStateIndex][amvs_q_index][bra_index] = value;

                                ++bra_index;
                            }
                        }
                    }
                }
                assert(max_bra_index == bra_index);
                
                ++ket_index;
            }
        }
    }

    const int max_ket_index = ket_contractScales.size()
        * ((angularMomentumQ +1) * (angularMomentumQ +2) * (angularMomentumQ +3) / 6);
    assert(ket_index == max_ket_index);

#ifdef DEBUG_HGP
    std::cerr << "transpose() end:: isCalcdERI_.size()=" << this->isCalcdERI_.size() << std::endl;
#endif // DEBUG_HGP
}


void DfEriEngine::calcERI(const int a_bar, const int b_bar,
                          const int a, const int b, const int p,
                          const int a_prime, const int b_prime, const int p_prime,
                          EriDataType* pERI)
{
    if ((a_bar < 0) || (b_bar < 0) ||
        (a < 0) || (b < 0) || (p < 0) ||
        (a_prime < 0) || (b_prime < 0) || (p_prime < 0)) {
        // do nothing because this state does not exist!
        return;
    }

#ifdef CHECK_MAX_COUNT
    this->max_a_bar_ = std::max(this->max_a_bar_, a_bar);
    this->max_b_bar_ = std::max(this->max_b_bar_, b_bar);
    this->max_a_ = std::max(this->max_a_, a);
    this->max_b_ = std::max(this->max_b_, b);
    this->max_p_ = std::max(this->max_p_, p);
    this->max_a_prime_ = std::max(this->max_a_prime_, a_prime);
    this->max_b_prime_ = std::max(this->max_b_prime_, b_prime);
    this->max_p_prime_ = std::max(this->max_p_prime_, p_prime);
#endif // CHECK_MAX_COUNT

    const ERI_State eriState(a_bar, b_bar, a, b, p, a_prime, b_prime, p_prime);
    
#ifdef DEBUG_HGP
    this->log_.debug(TlUtils::format("calcERI %s calcdERIs=%d",
                                     eriState.debugOut().c_str(),
                                     (int)this->isCalcdERI_.size()));
#endif // DEBUG_HGP

    if (this->isCalcdERI_[eriState] != true) {
        // recursive
        if (b_bar > 0) {
            // use eq.47
            this->calcERI(a_bar,   b_bar-1, a,   b,   p+1, a_prime,   b_prime,   p_prime, pERI);
            this->calcERI(a_bar+1, b_bar-1, a,   b,   p,   a_prime,   b_prime,   p_prime, pERI);
            this->calcERI(a_bar,   b_bar-1, a-1, b,   p,   a_prime,   b_prime,   p_prime, pERI);
            this->calcERI(a_bar,   b_bar-1, a,   b-1, p,   a_prime,   b_prime,   p_prime, pERI);
            this->ERI_EQ47(eriState, pERI);
        } else if (b > 0) {
            // use eq.44
            this->calcERI(a_bar,   b_bar,   a+1, b-1, p,   a_prime,   b_prime,   p_prime, pERI);
            this->calcERI(a_bar,   b_bar,   a,   b-1, p,   a_prime,   b_prime,   p_prime, pERI);
            this->calcERI(a_bar-1, b_bar,   a,   b-1, p,   a_prime,   b_prime,   p_prime, pERI);
            this->calcERI(a_bar,   b_bar-1, a,   b-1, p,   a_prime,   b_prime,   p_prime, pERI);
            this->ERI_EQ44(eriState, pERI);
        } else if (a_bar > 0) {
            // use eq.43
            this->calcERI(a_bar-1, b_bar,   a+1, b,   p,   a_prime+1, b_prime,   p_prime, pERI);
            this->calcERI(a_bar-1, b_bar,   a-1, b,   p,   a_prime,   b_prime,   p_prime, pERI);
            this->calcERI(a_bar-1, b_bar,   a,   b,   p-1, a_prime+1, b_prime,   p_prime, pERI);
            this->ERI_EQ43(eriState, pERI);
        } else if (a > 0) {
            // use eq.45
            this->calcERI(a_bar,   b_bar,   a-1, b,   p-1, a_prime,   b_prime,   p_prime  , pERI);
            this->calcERI(a_bar,   b_bar,   a-1, b,   p+1, a_prime,   b_prime,   p_prime+1, pERI);
            this->calcERI(a_bar,   b_bar,   a-1, b,   p,   a_prime,   b_prime+1, p_prime+1, pERI);
            this->calcERI(a_bar-1, b_bar,   a-1, b,   p,   a_prime,   b_prime+1, p_prime+1, pERI);
            this->calcERI(a_bar,   b_bar-1, a-1, b,   p+1, a_prime,   b_prime,   p_prime+1, pERI);
            this->ERI_EQ45(eriState, pERI);
        } else {
            // something wrong
            this->log_.critical("DfEriEngine::calcERI(): not calcd.");
            this->log_.critical(TlUtils::format(" [a~=%d, b~=%d, a=%d, b=%d, p=%d, a'=%d, b='%d, p='%d]",
                                                a_bar, b_bar, a, b, p, a_prime, b_prime, p_prime));
            abort();
        }
        
        this->isCalcdERI_[eriState] = true;
    }
}


int DfEriEngine::index(const TlAngularMomentumVector& a_bar, const TlAngularMomentumVector& b_bar,
                       const TlAngularMomentumVector& a, const TlAngularMomentumVector& b,
                       const TlAngularMomentumVector& p) const
{
    //const int am_a_bar = a_bar.angularMomentum();
    const int am_b_bar = b_bar.angularMomentum();
    const int am_a = a.angularMomentum();
    const int am_b = b.angularMomentum();
    const int am_p = p.angularMomentum();

    //const int a_bars = TlAngularMomentumVectorSet(am_a_bar).size();
    const int b_bars = TlAngularMomentumVectorSet(am_b_bar).size();
    const int as = TlAngularMomentumVectorSet(am_a).size();
    const int bs = TlAngularMomentumVectorSet(am_b).size();
    const int ps = TlAngularMomentumVectorSet(am_p).size();

    const int index = ((((a_bar.index()) * b_bars +b_bar.index()) * as +a.index()) * bs +b.index()) * ps +p.index();

    // assert(cs_index < (ERI_P_PRIME_MAX * ERI_P_PRIME_MAX * ERI_P_PRIME_MAX *
    //                    ERI_P_PRIME_MAX * ERI_P_PRIME_MAX * ERI_P_PRIME_MAX * ERI_NUM_OF_R_KINDS));
    return index;
}


// eq.44 (a -> b)
void DfEriEngine::ERI_EQ44(const ERI_State eriState, EriDataType* pERI)
{
#ifdef DEBUG_EQ44
    std::cerr << TlUtils::format("ERI_EQ44(a->b) entered: for %s",
                                 eriState.debugOut().c_str())
              << std::endl;
#endif // DEBUG_EQ44

    const int a_bar = eriState.a_bar;
    const int b_bar = eriState.b_bar;
    const int a = eriState.a;
    const int b = eriState.b;
    const int p = eriState.p;
    const int a_prime = eriState.a_prime;
    const int b_prime = eriState.b_prime;
    const int p_prime = eriState.p_prime;

    const TlAngularMomentumVectorSet AMVS_a_bar(a_bar);
    const TlAngularMomentumVectorSet AMVS_b_bar(b_bar);
    const TlAngularMomentumVectorSet AMVS_a(a);
    const TlAngularMomentumVectorSet AMVS_b(b);
    const TlAngularMomentumVectorSet AMVS_p(p);
    const int numOfAmv_a_bar = AMVS_a_bar.size();
    const int numOfAmv_b_bar = AMVS_b_bar.size();
    const int numOfAmv_a = AMVS_a.size();
    const int numOfAmv_b = AMVS_b.size();
    const int numOfAmv_p = AMVS_p.size();

    const int numOfElements = numOfAmv_a_bar * numOfAmv_b_bar * numOfAmv_a * numOfAmv_b * numOfAmv_p;
    const std::size_t stateIndex = eriState.index();
    (*pERI)[stateIndex].resize(numOfElements);

#ifdef CHECK_MAX_COUNT
    this->maxNumOfAMVs_ = std::max(this->maxNumOfAMVs_, numOfElements);
#endif // CHECK_MAX_COUNT

    for (int amv_b_index = 0; amv_b_index < numOfAmv_b; ++amv_b_index) {
        const TlAngularMomentumVector amv_b = AMVS_b.get(amv_b_index);
        const int i = this->initiativeRM(amv_b);
        const double AB_i = this->AB_[i];
        
        const TlAngularMomentumVector amv_b1 = amv_b - this->E1_[i];
        assert(amv_b1.isExist() == true);
        
        for (int amv_a_bar_index = 0; amv_a_bar_index < numOfAmv_a_bar; ++amv_a_bar_index) {
            const TlAngularMomentumVector amv_a_bar = AMVS_a_bar.get(amv_a_bar_index);
            const TlAngularMomentumVector amv_a_bar1 = amv_a_bar - this->E1_[i];
            
            for (int amv_b_bar_index = 0; amv_b_bar_index < numOfAmv_b_bar; ++amv_b_bar_index) {
                const TlAngularMomentumVector amv_b_bar = AMVS_b_bar.get(amv_b_bar_index);
                const TlAngularMomentumVector amv_b_bar1 = amv_b_bar - this->E1_[i];
            
                for (int amv_a_index = 0; amv_a_index < numOfAmv_a; ++amv_a_index) {
                    const TlAngularMomentumVector amv_a = AMVS_a.get(amv_a_index);
                    const TlAngularMomentumVector amv_a1p = amv_a + this->E1_[i];
                    
                    for (int amv_p_index = 0; amv_p_index < numOfAmv_p; ++amv_p_index) {
                        const TlAngularMomentumVector amv_p = AMVS_p.get(amv_p_index);

                        const std::size_t index = this->index(amv_a_bar, amv_b_bar,
                                                              amv_a, amv_b, amv_p);
                        //(*pERI)[stateIndex][index].resize(this->ERI_batch_);
                        (*pERI)[stateIndex][index].resize(ERI_MAX_BATCH);

                        // batch -----------------------------------------------
                        const ERI_State eriState1(a_bar, b_bar, a +1, b -1, p,
                                                  a_prime, b_prime, p_prime);
                        const std::size_t eriState1_index = eriState1.index();
                        const std::size_t index1 = this->index(amv_a_bar, amv_b_bar,
                                                               amv_a1p, amv_b1, amv_p);

                        const ERI_State eriState2(a_bar, b_bar, a, b -1, p,
                                                  a_prime, b_prime, p_prime);
                        const std::size_t eriState2_index = eriState2.index();
                        const std::size_t index2 = this->index(amv_a_bar, amv_b_bar,
                                                               amv_a, amv_b1, amv_p);

                        std::size_t eriState3_index = 0;
                        std::size_t index3 = 0;                        
                        if (amv_a_bar1.isExist() == true) {
                            const ERI_State eriState3(a_bar -1, b_bar, a, b -1, p,
                                                      a_prime, b_prime, p_prime);
                            eriState3_index = eriState3.index();
                            index3 = this->index(amv_a_bar1, amv_b_bar,
                                                 amv_a, amv_b1, amv_p);
                        }

                        std::size_t eriState4_index = 0;
                        std::size_t index4 = 0;
                        if (amv_b_bar1.isExist() == true) {
                            const ERI_State eriState4(a_bar, b_bar -1, a, b -1, p,
                                                      a_prime, b_prime, p_prime);
                            eriState4_index = eriState4.index();
                            index4 = this->index(amv_a_bar, amv_b_bar1,
                                                 amv_a, amv_b1, amv_p);
                        }
                        
                        for (int batch = 0; batch < this->ERI_batch_; ++batch) {

                            double answer = 0.0;
                            
                            // 1st term
                            {
                                const double value1 = (*pERI)[eriState1_index][index1][batch];
#ifdef DEBUG_EQ44
                                std::cerr << TlUtils::format("EQ44(1st): (%s %s, %s %s %s, %d %d %d| batch[%d] = % f",
                                                             amv_a_bar.debugOut().c_str(),
                                                             amv_b_bar.debugOut().c_str(),
                                                             amv_a1p.debugOut().c_str(),
                                                             amv_b1.debugOut().c_str(),
                                                             amv_p.debugOut().c_str(),
                                                             a_prime, b_prime, p_prime,
                                                             batch, value1)
                                          << std::endl;
#endif // DEBUG_EQ44
                                answer += value1;
                            }
                            
                            // 2nd term
                            {
                                const double value2 = (*pERI)[eriState2_index][index2][batch];
#ifdef DEBUG_EQ44
                                std::cerr << TlUtils::format("EQ44(2nd): (%s %s, %s %s %s, %d %d %d| batch[%d] = % f",
                                                             amv_a_bar.debugOut().c_str(),
                                                             amv_b_bar.debugOut().c_str(),
                                                             amv_a.debugOut().c_str(),
                                                             amv_b1.debugOut().c_str(),
                                                             amv_p.debugOut().c_str(),
                                                             a_prime, b_prime, p_prime,
                                                             batch, value2)
                                          << std::endl;
#endif // DEBUG_EQ44
                                // answer -= this->AB_[i] * value2;
                                answer -= AB_i * value2;
                            }
                            
                            // 3rd term
                            if (amv_a_bar1.isExist() == true) {
                                answer += amv_a_bar.get(i) * (*pERI)[eriState3_index][index3][batch];
                            }
                
                            // 4th term
                            if (amv_b_bar1.isExist() == true) {
                                answer -= amv_b_bar.get(i) * (*pERI)[eriState4_index][index4][batch];
                            }

                            (*pERI)[stateIndex][index][batch] = answer;
                        }
                    }
                }
            }
        }
    }
}


// eq.45 (p -> a)
void DfEriEngine::ERI_EQ45(const ERI_State eriState, EriDataType* pERI)
{
#ifdef DEBUG_EQ45
    std::cerr << TlUtils::format("ERI_EQ45(p->a) entered: for %s",
                                 eriState.debugOut().c_str())
              << std::endl;
#endif // DEBUG_EQ45
    
    const int a_bar = eriState.a_bar;
    const int b_bar = eriState.b_bar;
    const int a = eriState.a;
    const int b = eriState.b;
    const int p = eriState.p;
    const int a_prime = eriState.a_prime;
    const int b_prime = eriState.b_prime;
    const int p_prime = eriState.p_prime;

    const TlAngularMomentumVectorSet AMVS_a_bar(a_bar);
    const TlAngularMomentumVectorSet AMVS_b_bar(b_bar);
    const TlAngularMomentumVectorSet AMVS_a(a);
    const TlAngularMomentumVectorSet AMVS_b(b);
    const TlAngularMomentumVectorSet AMVS_p(p);
    const int numOfAmv_a_bar = AMVS_a_bar.size();
    const int numOfAmv_b_bar = AMVS_b_bar.size();
    const int numOfAmv_a = AMVS_a.size();
    const int numOfAmv_b = AMVS_b.size();
    const int numOfAmv_p = AMVS_p.size();

    const int numOfElements = numOfAmv_a_bar * numOfAmv_b_bar * numOfAmv_a * numOfAmv_b * numOfAmv_p;
    const std::size_t stateIndex = eriState.index();
    (*pERI)[stateIndex].resize(numOfElements);

#ifdef CHECK_MAX_COUNT
    this->maxNumOfAMVs_ = std::max(this->maxNumOfAMVs_, numOfElements);
#endif // CHECK_MAX_COUNT
    
    // p -> a
    for (int amv_a_index = 0; amv_a_index < numOfAmv_a; ++amv_a_index) {
        const TlAngularMomentumVector amv_a = AMVS_a.get(amv_a_index);
        const int i = this->initiativeRM(amv_a);
        const double AB_i = this->AB_[i];

        const TlAngularMomentumVector amv_a1 = amv_a - this->E1_[i];
        assert(amv_a1.isExist() == true);

        for (int amv_a_bar_index = 0; amv_a_bar_index < numOfAmv_a_bar; ++amv_a_bar_index) {
            const TlAngularMomentumVector amv_a_bar = AMVS_a_bar.get(amv_a_bar_index);
            const TlAngularMomentumVector amv_a_bar1 = amv_a_bar - this->E1_[i];

            for (int amv_b_bar_index = 0; amv_b_bar_index < numOfAmv_b_bar; ++amv_b_bar_index) {
                const TlAngularMomentumVector amv_b_bar = AMVS_b_bar.get(amv_b_bar_index);
                const TlAngularMomentumVector amv_b_bar1 = amv_b_bar - this->E1_[i];

                for (int amv_b_index = 0; amv_b_index < numOfAmv_b; ++amv_b_index) {
                    const TlAngularMomentumVector amv_b = AMVS_b.get(amv_b_index);

                    for (int amv_p_index = 0; amv_p_index < numOfAmv_p; ++amv_p_index) {
                        const TlAngularMomentumVector amv_p = AMVS_p.get(amv_p_index);
                        const TlAngularMomentumVector amv_p1 = amv_p - this->E1_[i];
                        const TlAngularMomentumVector amv_p1p = amv_p + this->E1_[i];
                            
                        const std::size_t index = this->index(amv_a_bar, amv_b_bar,
                                                              amv_a, amv_b, amv_p);
                        //(*pERI)[stateIndex][index].resize(this->ERI_batch_);
                        (*pERI)[stateIndex][index].resize(ERI_MAX_BATCH);

                        // betch -----------------------------------------------
                        std::size_t eriState1_index = 0;
                        std::size_t index1 = 0;
                        if (amv_p1.isExist() == true) {
                            const ERI_State eriState1(a_bar, b_bar, a -1, b, p -1,
                                                      a_prime, b_prime, p_prime);
                            eriState1_index = eriState1.index();
                            index1 = this->index(amv_a_bar, amv_b_bar,
                                                 amv_a1, amv_b, amv_p1);
                        }
                        const double p_i = amv_p.get(i);

                        const ERI_State eriState2(a_bar, b_bar, a -1, b, p +1,
                                                  a_prime, b_prime, p_prime +1);
                        const std::size_t eriState2_index = eriState2.index();
                        const std::size_t index2 = this->index(amv_a_bar, amv_b_bar,
                                                               amv_a1, amv_b, amv_p1p);

                        const ERI_State eriState3(a_bar, b_bar, a -1, b, p,
                                                  a_prime, b_prime +1, p_prime +1);
                        const std::size_t eriState3_index = eriState3.index();
                        const std::size_t index3 = this->index(amv_a_bar, amv_b_bar,
                                                               amv_a1, amv_b, amv_p);

                        std::size_t eriState4_index = 0;
                        std::size_t index4 = 0;
                        if (amv_a_bar1.isExist() == true) {
                            const ERI_State eriState4(a_bar -1, b_bar, a -1, b, p,
                                                      a_prime, b_prime +1, p_prime +1);
                            eriState4_index = eriState4.index();
                            index4 = this->index(amv_a_bar1, amv_b_bar,
                                                 amv_a1, amv_b, amv_p);
                        }

                        std::size_t eriState5_index = 0;
                        std::size_t index5 = 0;
                        if (amv_b_bar1.isExist() == true) {
                            const ERI_State eriState5(a_bar, b_bar -1, a -1, b, p,
                                                      a_prime, b_prime +1, p_prime +1);
                            eriState5_index = eriState5.index();
                            index5 = this->index(amv_a_bar, amv_b_bar1,
                                                 amv_a1, amv_b, amv_p);
                        }

                        const int batchSize = this->ERI_batch_;
                        for (int batch = 0; batch < batchSize; ++batch) {
                            double answer = 0.0;
                            
                            // 1st term
                            if (amv_p1.isExist() == true) {
                                const double value1 = (*pERI)[eriState1_index][index1][batch];
#ifdef DEBUG_EQ45
                                std::cerr << TlUtils::format("EQ45(1st): (%s %s, %s %s %s, %d %d %d| batch[%d] = % f",
                                                             amv_a_bar.debugOut().c_str(),
                                                             amv_b_bar.debugOut().c_str(),
                                                             amv_a1.debugOut().c_str(),
                                                             amv_b.debugOut().c_str(),
                                                             amv_p1.debugOut().c_str(),
                                                             a_prime, b_prime, p_prime,
                                                             batch, value1)
                                          << std::endl;
#endif // DEBUG_EQ45
                                answer += p_i * value1;
                            }
                            
                            // 2nd term
                            {
                                const double value2 = (*pERI)[eriState2_index][index2][batch];
#ifdef DEBUG_EQ45
                                std::cerr << TlUtils::format("EQ45(2nd): (%s %s, %s %s %s, %d %d %d| batch[%d] = % f",
                                                             amv_a_bar.debugOut().c_str(),
                                                             amv_b_bar.debugOut().c_str(),
                                                             amv_a1.debugOut().c_str(),
                                                             amv_b.debugOut().c_str(),
                                                             amv_p1p.debugOut().c_str(),
                                                             a_prime, b_prime, p_prime +1,
                                                             batch, value2)
                                          << std::endl;
#endif // DEBUG_EQ45
                                answer += value2;
                            }
                            
                            // 3rd term
                            {
                                const double value3 = (*pERI)[eriState3_index][index3][batch];
#ifdef DEBUG_EQ45
                                std::cerr << TlUtils::format("EQ45(3rd): (%s %s, %s %s %s, %d %d %d| batch[%d] = % f",
                                                             amv_a_bar.debugOut().c_str(),
                                                             amv_b_bar.debugOut().c_str(),
                                                             amv_a1.debugOut().c_str(),
                                                             amv_b.debugOut().c_str(),
                                                             amv_p.debugOut().c_str(),
                                                             a_prime, b_prime+1, p_prime+1,
                                                             batch, value3)
                                          << std::endl;
#endif // DEBUG_EQ45
                                // answer += this->AB_[i] * value3;
                                answer += AB_i * value3;
                            }
                            
                            // 4th term
                            if (amv_a_bar1.isExist() == true) {
                                answer -= amv_a_bar.get(i) * (*pERI)[eriState4_index][index4][batch];
                            }
                            
                            // 5th term
                            if (amv_b_bar1.isExist() == true) {
                                answer += amv_b_bar.get(i) * (*pERI)[eriState5_index][index5][batch];
                            }
                            
                            (*pERI)[stateIndex][index][batch] = answer;
                        }
                    }
                }
            }
        }
    }
}


// eq.46 (p -> a^)
void DfEriEngine::ERI_EQ46(const ERI_State eriState, EriDataType* pERI)
{
#ifdef DEBUG_EQ46
    std::cerr << TlUtils::format("ERI_EQ46(p->a^) entered: for %s",
                                 eriState.debugOut().c_str())
              << std::endl;
#endif // DEBUG_EQ46
    
    const int a_bar = eriState.a_bar;
    const int b_bar = eriState.b_bar;
    const int a = eriState.a;
    const int b = eriState.b;
    const int p = eriState.p;
    const int a_prime = eriState.a_prime;
    const int b_prime = eriState.b_prime;
    const int p_prime = eriState.p_prime;

    const TlAngularMomentumVectorSet AMVS_a_bar(a_bar);
    const TlAngularMomentumVectorSet AMVS_b_bar(b_bar);
    const TlAngularMomentumVectorSet AMVS_a(a);
    const TlAngularMomentumVectorSet AMVS_b(b);
    const TlAngularMomentumVectorSet AMVS_p(p);
    const int numOfAmv_a_bar = AMVS_a_bar.size();
    const int numOfAmv_b_bar = AMVS_b_bar.size();
    const int numOfAmv_a = AMVS_a.size();
    const int numOfAmv_b = AMVS_b.size();
    const int numOfAmv_p = AMVS_p.size();

    const int numOfElements = numOfAmv_a_bar * numOfAmv_b_bar * numOfAmv_a * numOfAmv_b * numOfAmv_p;
    (*pERI)[eriState.index()].resize(numOfElements);

#ifdef CHECK_MAX_COUNT
    this->maxNumOfAMVs_ = std::max(this->maxNumOfAMVs_, numOfElements);
#endif // CHECK_MAX_COUNT

    // p -> a^
    for (int amv_a_bar_index = 0; amv_a_bar_index < numOfAmv_a_bar; ++amv_a_bar_index) {
        const TlAngularMomentumVector amv_a_bar = AMVS_a_bar.get(amv_a_bar_index);
        const int i = this->initiativeRM(amv_a_bar);

        for (int amv_b_bar_index = 0; amv_b_bar_index < numOfAmv_b_bar; ++amv_b_bar_index) {
            const TlAngularMomentumVector amv_b_bar = AMVS_b_bar.get(amv_b_bar_index);

            for (int amv_a_index = 0; amv_a_index < numOfAmv_a; ++amv_a_index) {
                const TlAngularMomentumVector amv_a = AMVS_a.get(amv_a_index);

                for (int amv_b_index = 0; amv_b_index < numOfAmv_b; ++amv_b_index) {
                    const TlAngularMomentumVector amv_b = AMVS_b.get(amv_b_index);

                    for (int amv_p_index = 0; amv_p_index < numOfAmv_p; ++amv_p_index) {
                        const TlAngularMomentumVector amv_p = AMVS_p.get(amv_p_index);

                        const std::size_t index = this->index(amv_a_bar, amv_b_bar,
                                                              amv_a, amv_b, amv_p);
                        //(*pERI)[eriState.index()][index].resize(this->ERI_batch_);
                        (*pERI)[eriState.index()][index].resize(ERI_MAX_BATCH);

                        // betch -----------------------------------------------
                        for (int batch = 0; batch < this->ERI_batch_; ++batch) {
                            double answer = 0.0;
                            
                            const TlAngularMomentumVector amv_a_bar1 = amv_a_bar - this->E1_[i];
                            assert(amv_a_bar1.isExist() == true);
                            
                            // 1st term
                            {
                                const TlAngularMomentumVector amv_p1p = amv_p + this->E1_[i];
                                const ERI_State eriState1(a_bar -1, b_bar, a, b, p +1,
                                                          a_prime +1, b_prime, p_prime +1);
                                const std::size_t index1 = this->index(amv_a_bar1, amv_b_bar,
                                                                       amv_a, amv_b, amv_p1p);
                                //const double value1 = (*pERI)[eriState1][index1][batch];
                                const double value1 = (*pERI)[eriState1.index()][index1][batch];
#ifdef DEBUG_EQ46
                                std::cerr << TlUtils::format("EQ46(1st): (%s %s, %s %s %s, %d %d %d| batch[%d] = % f",
                                                             amv_a_bar1.debugOut().c_str(),
                                                             amv_b_bar.debugOut().c_str(),
                                                             amv_a.debugOut().c_str(),
                                                             amv_b.debugOut().c_str(),
                                                             amv_p1p.debugOut().c_str(),
                                                             a_prime +1, b_prime, p_prime +1,
                                                             batch, value1)
                                          << std::endl;
#endif // DEBUG_EQ46
                                answer += value1;
                            }
                            
                            // 2nd term
                            const TlAngularMomentumVector amv_a1 = amv_a - this->E1_[i];
                            if (amv_a1.isExist() == true) {
                                const ERI_State eriState2(a_bar -1, b_bar, a -1, b, p,
                                                          a_prime, b_prime, p_prime);
                                const std::size_t index2 = this->index(amv_a_bar1, amv_b_bar,
                                                                       amv_a1, amv_b, amv_p);
                                //const double value2 = (*pERI)[eriState2][index2][batch];
                                const double value2 = (*pERI)[eriState2.index()][index2][batch];
#ifdef DEBUG_EQ46
                                std::cerr << TlUtils::format("EQ46(2nd): (%s %s, %s %s %s, %d %d %d| batch[%d] = % f",
                                                             amv_a_bar1.debugOut().c_str(),
                                                             amv_b_bar.debugOut().c_str(),
                                                             amv_a1.debugOut().c_str(),
                                                             amv_b.debugOut().c_str(),
                                                             amv_p.debugOut().c_str(),
                                                             a_prime, b_prime, p_prime,
                                                             batch, value2)
                                          << std::endl;
#endif // DEBUG_EQ46
                                answer -= amv_a.get(i) * value2;
                            }
                            
                            // 3rd term
                            {
                                const ERI_State eriState3(a_bar -1, b_bar, a, b, p,
                                                          a_prime +1, b_prime +1, p_prime +1);
                                const std::size_t index3 = this->index(amv_a_bar1, amv_b_bar,
                                                                       amv_a, amv_b, amv_p);
                                //const double value3 = (*pERI)[eriState3][index3][batch];
                                const double value3 = (*pERI)[eriState3.index()][index3][batch];
#ifdef DEBUG_EQ46
                                std::cerr << TlUtils::format("EQ46(3rd): (%s %s, %s %s %s, %d %d %d| batch[%d] = % f",
                                                             amv_a_bar1.debugOut().c_str(),
                                                             amv_b_bar.debugOut().c_str(),
                                                             amv_a.debugOut().c_str(),
                                                             amv_b.debugOut().c_str(),
                                                             amv_p.debugOut().c_str(),
                                                             a_prime +1, b_prime +1, p_prime +1,
                                                             batch, value3)
                                          << std::endl;
#endif // DEBUG_EQ46
                                answer += this->AB_[i] * value3;
                            }
                            
                            // 4th term
                            const TlAngularMomentumVector amv_a_bar2 = amv_a_bar - this->E2_[i];
                            if (amv_a_bar2.isExist() == true) {
                                const ERI_State eriState4(a_bar -2, b_bar, a, b, p,
                                                          a_prime +1, b_prime +1, p_prime +1);
                                const std::size_t index4 = this->index(amv_a_bar2, amv_b_bar,
                                                                       amv_a, amv_b, amv_p);
                                //answer -= (amv_a_bar.get(i) -1.0) * (*pERI)[eriState4][index4][batch];
                                answer -= (amv_a_bar.get(i) -1.0) * (*pERI)[eriState4.index()][index4][batch];
                            }
                            
                            // 5th term
                            const TlAngularMomentumVector amv_b_bar1 = amv_b_bar - this->E1_[i];
                            if (amv_b_bar1.isExist() == true) {
                                const ERI_State eriState5(a_bar -1, b_bar -1, a, b, p,
                                                          a_prime +1, b_prime +1, p_prime +1);
                                const std::size_t index5 = this->index(amv_a_bar1, amv_b_bar1,
                                                                       amv_a, amv_b, amv_p);
                                //answer += amv_b_bar.get(i) * (*pERI)[eriState5][index5][batch];
                                answer += amv_b_bar.get(i) * (*pERI)[eriState5.index()][index5][batch];
                            }
                            
                            //(*pERI)[eriState][index][batch] = answer;
                            (*pERI)[eriState.index()][index][batch] = answer;
                        }
                    }
                }
            }
        }
    }
}


// eq.43 (a -> a^)
void DfEriEngine::ERI_EQ43(const ERI_State eriState, EriDataType* pERI)
{
#ifdef DEBUG_EQ43
    std::cerr << TlUtils::format("ERI_EQ43(a->a^) entered: for %s",
                                 eriState.debugOut().c_str())
              << std::endl;
#endif // DEBUG_EQ43
    
    const int a_bar = eriState.a_bar;
    const int b_bar = eriState.b_bar;
    const int a = eriState.a;
    const int b = eriState.b;
    const int p = eriState.p;
    const int a_prime = eriState.a_prime;
    const int b_prime = eriState.b_prime;
    const int p_prime = eriState.p_prime;

    const TlAngularMomentumVectorSet AMVS_a_bar(a_bar);
    const TlAngularMomentumVectorSet AMVS_b_bar(b_bar);
    const TlAngularMomentumVectorSet AMVS_a(a);
    const TlAngularMomentumVectorSet AMVS_b(b);
    const TlAngularMomentumVectorSet AMVS_p(p);
    const int numOfAmv_a_bar = AMVS_a_bar.size();
    const int numOfAmv_b_bar = AMVS_b_bar.size();
    const int numOfAmv_a = AMVS_a.size();
    const int numOfAmv_b = AMVS_b.size();
    const int numOfAmv_p = AMVS_p.size();

    const int numOfElements = numOfAmv_a_bar * numOfAmv_b_bar * numOfAmv_a * numOfAmv_b * numOfAmv_p;
    (*pERI)[eriState.index()].resize(numOfElements);
    
#ifdef CHECK_MAX_COUNT
    this->maxNumOfAMVs_ = std::max(this->maxNumOfAMVs_, numOfElements);
#endif // CHECK_MAX_COUNT

    // a -> a^
    for (int amv_a_bar_index = 0; amv_a_bar_index < numOfAmv_a_bar; ++amv_a_bar_index) {
        const TlAngularMomentumVector amv_a_bar = AMVS_a_bar.get(amv_a_bar_index);
        const int i = this->initiativeRM(amv_a_bar);

        for (int amv_b_bar_index = 0; amv_b_bar_index < numOfAmv_b_bar; ++amv_b_bar_index) {
            const TlAngularMomentumVector amv_b_bar = AMVS_b_bar.get(amv_b_bar_index);

            for (int amv_a_index = 0; amv_a_index < numOfAmv_a; ++amv_a_index) {
                const TlAngularMomentumVector amv_a = AMVS_a.get(amv_a_index);

                for (int amv_b_index = 0; amv_b_index < numOfAmv_b; ++amv_b_index) {
                    const TlAngularMomentumVector amv_b = AMVS_b.get(amv_b_index);

                    for (int amv_p_index = 0; amv_p_index < numOfAmv_p; ++amv_p_index) {
                        const TlAngularMomentumVector amv_p = AMVS_p.get(amv_p_index);

                        const std::size_t index = this->index(amv_a_bar, amv_b_bar,
                                                              amv_a, amv_b, amv_p);
                        //(*pERI)[eriState.index()][index].resize(this->ERI_batch_);
                        (*pERI)[eriState.index()][index].resize(ERI_MAX_BATCH);

                        // betch -----------------------------------------------
                        for (int batch = 0; batch < this->ERI_batch_; ++batch) {
                            double answer = 0.0;
                            
                            const TlAngularMomentumVector amv_a_bar1 = amv_a_bar - this->E1_[i];
                            assert(amv_a_bar1.isExist() == true);
                            
                            // 1st term
                            {
                                const TlAngularMomentumVector amv_a1p = amv_a + this->E1_[i];
                                const ERI_State eriState1(a_bar -1, b_bar, a +1, b, p,
                                                          a_prime +1, b_prime, p_prime);
                                const std::size_t index1 = this->index(amv_a_bar1, amv_b_bar,
                                                                       amv_a1p, amv_b, amv_p);
                                //const double value1 = (*pERI)[eriState1][index1][batch];
                                const double value1 = (*pERI)[eriState1.index()][index1][batch];
#ifdef DEBUG_EQ43
                                std::cerr << TlUtils::format("EQ43(1st): (%s %s, %s %s %s, %d %d %d| batch[%d] = % f",
                                                             amv_a_bar1.debugOut().c_str(),
                                                             amv_b_bar.debugOut().c_str(),
                                                             amv_a1p.debugOut().c_str(),
                                                             amv_b.debugOut().c_str(),
                                                             amv_p.debugOut().c_str(),
                                                             a_prime +1, b_prime, p_prime,
                                                             batch, value1)
                                          << std::endl;
#endif // DEBUG_EQ43
                                answer += value1;
                            }
                            
                            // 2nd term
                            const TlAngularMomentumVector amv_a1 = amv_a - this->E1_[i];
                            if (amv_a1.isExist() == true) {
                                const ERI_State eriState2(a_bar -1, b_bar, a -1, b, p,
                                                          a_prime, b_prime, p_prime);
                                const std::size_t index2 = this->index(amv_a_bar1, amv_b_bar,
                                                                       amv_a1, amv_b, amv_p);
                                //const double value2 = (*pERI)[eriState2][index2][batch];
                                const double value2 = (*pERI)[eriState2.index()][index2][batch];
#ifdef DEBUG_EQ43
                                std::cerr << TlUtils::format("EQ43(2nd): (%s %s, %s %s %s, %d %d %d| batch[%d] = % f",
                                                             amv_a_bar1.debugOut().c_str(),
                                                             amv_b_bar.debugOut().c_str(),
                                                             amv_a1.debugOut().c_str(),
                                                             amv_b.debugOut().c_str(),
                                                             amv_p.debugOut().c_str(),
                                                             a_prime, b_prime, p_prime,
                                                             batch, value2)
                                          << std::endl;
#endif // DEBUG_EQ43
                                answer -= amv_a.get(i) * value2;
                            }
                            
                            // 3rd term
                            const TlAngularMomentumVector amv_p1 = amv_p - this->E1_[i];
                            if (amv_p1.isExist() == true) {
                                const ERI_State eriState3(a_bar -1, b_bar, a, b, p -1,
                                                          a_prime +1, b_prime, p_prime);
                                const std::size_t index3 = this->index(amv_a_bar1, amv_b_bar,
                                                                       amv_a, amv_b, amv_p1);
                                //const double value3 = (*pERI)[eriState3][index3][batch];
                                const double value3 = (*pERI)[eriState3.index()][index3][batch];
#ifdef DEBUG_EQ43
                                std::cerr << TlUtils::format("EQ43(3rd): (%s %s, %s %s %s, %d %d %d| batch[%d] = % f",
                                                             amv_a_bar1.debugOut().c_str(),
                                                             amv_b_bar.debugOut().c_str(),
                                                             amv_a.debugOut().c_str(),
                                                             amv_b.debugOut().c_str(),
                                                             amv_p1.debugOut().c_str(),
                                                             a_prime, b_prime, p_prime,
                                                             batch, value3)
                                          << std::endl;
#endif // DEBUG_EQ43
                                answer -= amv_p.get(i) * value3;
                            }
                            
                            //(*pERI)[eriState][index][batch] = answer;
                            (*pERI)[eriState.index()][index][batch] = answer;
                        }
                    }
                }
            }
        }
    }
}


// eq.47 (a^ -> b^)
void DfEriEngine::ERI_EQ47(const ERI_State eriState, EriDataType* pERI)
{
#ifdef DEBUG_EQ47
    std::cerr << TlUtils::format("ERI_EQ47(a^->b^) entered: for %s",
                                 eriState.debugOut().c_str())
              << std::endl;
#endif // DEBUG_EQ47
    
    const int a_bar = eriState.a_bar;
    const int b_bar = eriState.b_bar;
    const int a = eriState.a;
    const int b = eriState.b;
    const int p = eriState.p;
    const int a_prime = eriState.a_prime;
    const int b_prime = eriState.b_prime;
    const int p_prime = eriState.p_prime;

    const TlAngularMomentumVectorSet AMVS_a_bar(a_bar);
    const TlAngularMomentumVectorSet AMVS_b_bar(b_bar);
    const TlAngularMomentumVectorSet AMVS_a(a);
    const TlAngularMomentumVectorSet AMVS_b(b);
    const TlAngularMomentumVectorSet AMVS_p(p);
    const int numOfAmv_a_bar = AMVS_a_bar.size();
    const int numOfAmv_b_bar = AMVS_b_bar.size();
    const int numOfAmv_a = AMVS_a.size();
    const int numOfAmv_b = AMVS_b.size();
    const int numOfAmv_p = AMVS_p.size();

    const int numOfElements = numOfAmv_a_bar * numOfAmv_b_bar * numOfAmv_a * numOfAmv_b * numOfAmv_p;
    (*pERI)[eriState.index()].resize(numOfElements);

#ifdef CHECK_MAX_COUNT
    this->maxNumOfAMVs_ = std::max(this->maxNumOfAMVs_, numOfElements);
#endif // CHECK_MAX_COUNT
    
    // a^ -> b^
    for (int amv_a_bar_index = 0; amv_a_bar_index < numOfAmv_a_bar; ++amv_a_bar_index) {
        const TlAngularMomentumVector amv_a_bar = AMVS_a_bar.get(amv_a_bar_index);

        for (int amv_b_bar_index = 0; amv_b_bar_index < numOfAmv_b_bar; ++amv_b_bar_index) {
            const TlAngularMomentumVector amv_b_bar = AMVS_b_bar.get(amv_b_bar_index);
            const int i = this->initiativeRM(amv_b_bar);

            for (int amv_a_index = 0; amv_a_index < numOfAmv_a; ++amv_a_index) {
                const TlAngularMomentumVector amv_a = AMVS_a.get(amv_a_index);

                for (int amv_b_index = 0; amv_b_index < numOfAmv_b; ++amv_b_index) {
                    const TlAngularMomentumVector amv_b = AMVS_b.get(amv_b_index);

                    for (int amv_p_index = 0; amv_p_index < numOfAmv_p; ++amv_p_index) {
                        const TlAngularMomentumVector amv_p = AMVS_p.get(amv_p_index);

                        const std::size_t index = this->index(amv_a_bar, amv_b_bar,
                                                              amv_a, amv_b, amv_p);
                        //(*pERI)[eriState.index()][index].resize(this->ERI_batch_);
                        (*pERI)[eriState.index()][index].resize(ERI_MAX_BATCH);

                        // betch -----------------------------------------------
                        for (int batch = 0; batch < this->ERI_batch_; ++batch) {
                            double answer = 0.0;
                            
                            const TlAngularMomentumVector amv_b_bar1 = amv_b_bar - this->E1_[i];
                            assert(amv_b_bar1.isExist() == true);
                            
                            // 1st term
                            {
                                const TlAngularMomentumVector amv_p1p = amv_p + this->E1_[i];
                                const ERI_State eriState1(a_bar, b_bar -1, a, b, p +1,
                                                          a_prime, b_prime, p_prime);
                                const std::size_t index1 = this->index(amv_a_bar, amv_b_bar1,
                                                                       amv_a, amv_b, amv_p1p);
                                const double value1 = (*pERI)[eriState1.index()][index1][batch];
#ifdef DEBUG_EQ47
                                std::cerr << TlUtils::format("EQ47(1st): (%s %s, %s %s %s, %d %d %d| batch[%d] = % f",
                                                             amv_a_bar.debugOut().c_str(),
                                                             amv_b_bar1.debugOut().c_str(),
                                                             amv_a.debugOut().c_str(),
                                                             amv_b.debugOut().c_str(),
                                                             amv_p1p.debugOut().c_str(),
                                                             a_prime, b_prime, p_prime,
                                                             batch, value1)
                                          << std::endl;
#endif // DEBUG_EQ47
                                answer += value1;
                            }
                            
                            // 2nd term
                            const TlAngularMomentumVector amv_a_bar1p = amv_a_bar + this->E1_[i];
                            {
                                const ERI_State eriState2(a_bar +1, b_bar -1, a, b, p,
                                                          a_prime, b_prime, p_prime);
                                const std::size_t index2 = this->index(amv_a_bar1p, amv_b_bar1,
                                                                       amv_a, amv_b, amv_p);
                                const double value2 = (*pERI)[eriState2.index()][index2][batch];
#ifdef DEBUG_EQ47
                                std::cerr << TlUtils::format("EQ47(2nd): (%s %s, %s %s %s, %d %d %d| batch[%d] = % f",
                                                             amv_a_bar1p.debugOut().c_str(),
                                                             amv_b_bar1.debugOut().c_str(),
                                                             amv_a.debugOut().c_str(),
                                                             amv_b.debugOut().c_str(),
                                                             amv_p.debugOut().c_str(),
                                                             a_prime, b_prime, p_prime,
                                                             batch, value2)
                                          << std::endl;
#endif // DEBUG_EQ47
                                answer -= value2;
                            }
                            
                            // 3rd term
                            const TlAngularMomentumVector amv_a1 = amv_a - this->E1_[i];
                            if (amv_a1.isExist() == true) {
                                const ERI_State eriState3(a_bar, b_bar -1, a -1, b, p,
                                                          a_prime, b_prime, p_prime);
                                const std::size_t index3 = this->index(amv_a_bar, amv_b_bar1,
                                                                       amv_a1, amv_b, amv_p);
                                const double value3 = (*pERI)[eriState3.index()][index3][batch];
#ifdef DEBUG_EQ47
                                std::cerr << TlUtils::format("EQ47(3rd): (%s %s, %s %s %s, %d %d %d| batch[%d] = % f",
                                                             amv_a_bar.debugOut().c_str(),
                                                             amv_b_bar1.debugOut().c_str(),
                                                             amv_a1.debugOut().c_str(),
                                                             amv_b.debugOut().c_str(),
                                                             amv_p.debugOut().c_str(),
                                                             a_prime, b_prime, p_prime,
                                                             batch, value3)
                                          << std::endl;
#endif // DEBUG_EQ46
                                answer -= amv_a.get(i) * value3;
                            }
                            
                            // 4th term
                            const TlAngularMomentumVector amv_b1 = amv_b - this->E1_[i];
                            if (amv_b1.isExist() == true) {
                                const ERI_State eriState4(a_bar, b_bar -1, a, b -1, p,
                                                          a_prime, b_prime, p_prime);
                                const std::size_t index4 = this->index(amv_a_bar, amv_b_bar1,
                                                                       amv_a, amv_b1, amv_p);
                                //answer -= amv_b.get(i) * (*pERI)[eriState4][index4][batch];
                                answer -= amv_b.get(i) * (*pERI)[eriState4.index()][index4][batch];
                            }
                            
                            //(*pERI)[eriState][index][batch] = answer;
                            (*pERI)[eriState.index()][index][batch] = answer;
                        }
                    }
                }
            }
        }
    }
}

