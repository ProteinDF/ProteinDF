#include <ios>
#include <iostream>
#include <fstream>
#include <cmath>
#include <utility>

#include "DfTwoElectronIntegral.h"
#include "TlPosition.h"
#include "TlFmt.h"
#include "TlUtils.h"
#include "TlTime.h"

#include "DfTEI.h"
#include "DfEriEngine.h"

////////////////////////////////////////////////////////////////////////
// static member variables
//
const double DfTwoElectronIntegral::CK          = sqrt(2.0 * sqrt(M_PI)) * M_PI;
const int DfTwoElectronIntegral::MAX_SHELL_TYPE = 3;
const double DfTwoElectronIntegral::CONTRIBUTE_COEF = 2.0 * std::pow(M_PI, 2.5);


////////////////////////////////////////////////////////////////////////
// construct and destruct
//
DfTwoElectronIntegral::DfTwoElectronIntegral(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam), m_TlOrbitalInfo((*pPdfParam)["coordinates"],
                                           (*pPdfParam)["basis_sets"])
{
    const TlSerializeData& pdfParam = *pPdfParam;
    // initialize
    this->m_dInputCutoffThreshold = pdfParam["cut-value"].getDouble();
    this->m_dCutoffThreshold = this->m_dInputCutoffThreshold;

    this->m_dPrimitiveIntegralsCutoffThreshold = 1.0E-16; //this->m_dInputCutoffThreshold * 1.0E-2;
    if (pdfParam["TEI-primitive-cutoff"].getStr() != "") {
        this->m_dPrimitiveIntegralsCutoffThreshold = pdfParam["TEI-primitive-cutoff"].getDouble(); //1.0E-15;
    }

    this->isDensityCutoff_ = true;
    if (TlUtils::toUpper(pdfParam["TEI-density-cutoff"].getStr()) == "NO") {
        this->isDensityCutoff_ = false;
    }

    this->lengthScaleParameter_ = 1.0;
    if (pdfParam["TEI-length-scale-parameter"].getStr() != "") {
        this->lengthScaleParameter_ = pdfParam["TEI-length-scale-parameter"].getDouble();
    }

    // for debug
    this->isOutputTEI_ = (TlUtils::toUpper(pdfParam["output_TEI"].getStr()) == "YES") ? true : false;
    this->outputTEI_path_ = "TEI.txt";
//     this->isOutputHGP_ = (TlUtils::toUpper(pdfParam["output_HGP"].getStr()) == "YES") ? true : false;

    //this->isUseNewEngine_ = (*pPdfParam)["new_engine"].getBoolean();
    
    //
    this->JShellListSize_ = 0;
    this->pJShellList_ = new index_type[this->m_nNumOfAOs];
    this->LShellListSize_ = 0;
    this->pLShellList_ = new index_type[this->m_nNumOfAOs];
}


DfTwoElectronIntegral::~DfTwoElectronIntegral()
{
    if (this->pJShellList_ != NULL) {
        delete[] this->pJShellList_;
        this->pJShellList_ = NULL;
    }
    if (this->pLShellList_ != NULL) {
        delete[] this->pLShellList_;
        this->pLShellList_ = NULL;
    }
}

void DfTwoElectronIntegral::prepare_ERI()
{
    // initialize
    this->m_ShellList.clear();
    this->m_DistributionShellList.clear();
    this->m_IKShellPairList.clear();

    // define threshold
    this->m_dCutoffThreshold = 1.0E-16;    

    // list up shell-pair
    this->makeShellList();
    this->screenDistributions();
    this->getShellList_density_nocut();
}


void DfTwoElectronIntegral::getContractKMatrixByIntegralDriven(const TlSymmetricMatrix& P,
                                                               TlSymmetricMatrix* pK)
{
    assert(pK != NULL);

    // initialize
    this->m_ShellList.clear();
    this->m_DistributionShellList.clear();
    this->m_IKShellPairList.clear();

    // define threshold
    {
        const double maxDensityElementValue = P.getMaxAbsoluteElement();
        this->m_dCutoffThreshold = this->m_dInputCutoffThreshold;
        if ((0.0 < maxDensityElementValue) && (maxDensityElementValue < 1.0)) {
            this->m_dCutoffThreshold /= maxDensityElementValue;
            this->logger(TlUtils::format(" cutoff value = %e\n", this->m_dCutoffThreshold));
        }
    }

    // list up shell-pair
    this->makeShellList();
    this->screenDistributions();
    if (this->isDensityCutoff_ == true) {
        this->getShellList_density();
    } else {
        this->getShellList_density_nocut();
    }

    this->getContractKMatrixByIntegralDrivenCore(P, pK);
}


void DfTwoElectronIntegral::getContractKMatrixByRTmethod(const TlSymmetricMatrix& P,
                                                         TlSymmetricMatrix* pK)
{

    assert(pK != NULL);

    // initialize
    this->m_ShellList.clear();
    this->m_DistributionShellList.clear();
    this->m_IKShellPairList.clear();

    // define threshold
    {
        const double maxDensityElementValue = P.getMaxAbsoluteElement();
        this->m_dCutoffThreshold = this->m_dInputCutoffThreshold;
        if ((0.0 < maxDensityElementValue) && (maxDensityElementValue < 1.0)) {
            this->m_dCutoffThreshold /= maxDensityElementValue;
            this->logger(TlUtils::format(" cutoff value = %e\n", this->m_dCutoffThreshold));
        }
    }

    // list up shell-pair
    this->makeShellList();
    this->screenDistributions();
    if (this->isDensityCutoff_ == true) {
        this->getShellList_density();
    } else {
        this->getShellList_density_nocut();
    }

    this->getContractKMatrixByRTmethodCore(P, pK);
}


// ---------------------------------------------------------------------
// integral driven
// ---------------------------------------------------------------------
void DfTwoElectronIntegral::getContractKMatrixByIntegralDrivenCore(const TlSymmetricMatrix& P,
                                                                   TlSymmetricMatrix* pK)
{
    assert(pK != NULL);

    const int nMaxIShell = this->m_TlOrbitalInfo.getNumOfOrbitals();
    pK->resize(nMaxIShell);

    // I --------------
    for (int iShellType = MAX_SHELL_TYPE -1; iShellType >= 0; --iShellType) {
        const std::vector<int> IShellList = this->m_ShellList[iShellType];

        // K --------------
        for (int kShellType = MAX_SHELL_TYPE -1; kShellType >= 0; --kShellType) {

            const int nIKShellPairIndex = 4*iShellType +kShellType;
            const std::vector<IKShellPair> aIKShellPairList = this->getLocalIKShellPairList(nIKShellPairIndex);
            std::vector<IKShellPair>::const_iterator ikShellPairEnd = aIKShellPairList.end();

            // J --------------
            for (int jShellType = MAX_SHELL_TYPE -1; jShellType >= 0; --jShellType) {

                // L --------------
                for (int lShellType = MAX_SHELL_TYPE -1; lShellType >= 0; --lShellType) {

                    const unsigned int eriType = DfTEI::getEriType(iShellType, jShellType, kShellType, lShellType);

                    for (std::vector<IKShellPair>::const_iterator ikShellPair = aIKShellPairList.begin();
                            ikShellPair != ikShellPairEnd; ++ikShellPair) {
                        const int iShell = ikShellPair->nIShell;
                        const int kShell = ikShellPair->nKShell;

                        this->IntegralDrivenCore(eriType, iShell, kShell, P, pK);
                    }
                }
            }
        }
    }

    this->finalize(*pK);
}


void DfTwoElectronIntegral::IntegralDrivenCore(const unsigned int eriType,
                                               const int iShell, const int kShell,
                                               const TlSymmetricMatrix& P,
                                               TlSymmetricMatrix* pK)
{
    assert(pK != NULL);

    // shell type
    const int iShellType = int((eriType >> 6) & 3);
    const int jShellType = int((eriType >> 4) & 3);
    const int kShellType = int((eriType >> 2) & 3);
    const int lShellType = int(eriType        & 3);

    // steps
    const int iSteps = 2*iShellType +1;
    const int jSteps = 2*jShellType +1;
    const int kSteps = 2*kShellType +1;
    const int lSteps = 2*lShellType +1;

    // j-Shell
    // iShellに対するjShellを受け取る。
//     const std::vector<int> JShellList
//     = this->getShellList_distribute(this->getShellList(this->m_ShellList[jShellType], iShell +1),
//                                     iShell, jShellType);
    // 行列の次元数はINT_MAXが最大なので、
    // 受け取る要素数(indexJ_max)はINT_MAX以下になる。
//     const int indexJ_max = JShellList.size();
    const int indexJ_max = this->getShellList_distribute(this->getShellList(this->m_ShellList[jShellType], iShell +1),
                                                         iShell, jShellType,
                                                         this->pJShellList_);

#pragma omp parallel
    {
        DfTEI dfTEI(this->m_dPrimitiveIntegralsCutoffThreshold);
        DfEriEngine dfEriEngine;
        
#pragma omp for schedule(runtime)
        for (int indexJ = 0; indexJ < indexJ_max; ++indexJ) {
            //const int jShell = JShellList[indexJ];
            const index_type jShell = this->pJShellList_[indexJ];
            const DfTEI::ShellPair IJ = this->getShellPair(iShell, jShell);
            const double dIJIJ = IJ.dSchwartzValue;
            //const double dP_jk = this->getMaxDensityMatrixElement(jShell, jSteps, kShell, kSteps, P);

            // l-Shell
            const int maxLShell = (kShell == iShell) ? (jShell) : (kShell);
//             const std::vector<int> LShellList =
//                 this->getShellList_distribute(this->getShellList(this->m_ShellList[lShellType], maxLShell +1),
//                                               kShell, lShellType);
//             const std::size_t indexL_max = LShellList.size();
            const index_type indexL_max = this->getShellList_distribute(this->getShellList(this->m_ShellList[lShellType], maxLShell +1),
                                                                        kShell, lShellType,
                                                                        this->pLShellList_);
            for (index_type indexL = 0; indexL < indexL_max; ++indexL) {
                //const int lShell = LShellList[indexL];
                const index_type lShell = this->pLShellList_[indexL];
                const DfTEI::ShellPair KL = this->getShellPair(kShell, lShell);
                const double dKLKL = KL.dSchwartzValue;
                //const double dP_il = this->getMaxDensityMatrixElement(iShell, iSteps, lShell, lSteps, P);
                //const double dP_jl = this->getMaxDensityMatrixElement(jShell, jSteps, lShell, lSteps, P);
                //const double dP_ik = this->getMaxDensityMatrixElement(iShell, iSteps, kShell, kSteps, P);

                // shwartz cutoff
                //double dMaxP = std::max(dP_jk, dP_il);
                //dMaxP = std::max(dMaxP, dP_jl);
                //dMaxP = std::max(dMaxP, dP_ik);
                //if ((dMaxP * dIJIJ * dKLKL) < this->m_dCutoffThreshold) {
                if (dIJIJ * dKLKL < this->m_dCutoffThreshold) {
                    continue;
                }

                // calc
//                 if (this->isUseNewEngine_ != true) {
//                     std::cerr << TlUtils::format("ERI calc: (%d %d | %d %d)",
//                                                  iShell, jShell, kShell, lShell)
//                               << std::endl;
                    dfTEI.calc(eriType, IJ, KL);
                    
                    // store K matrix
#pragma omp critical (storeK_ID)
                    {
                        this->storeKmatrixByIntegralDriven(iShell, iSteps, jShell, jSteps,
                                                           kShell, kSteps, lShell, lSteps, P, dfTEI.ERI, pK);
                    }
//                 } else {
//                     const DfEriEngine::Query queryPQ(0, 0, iShellType, jShellType);
//                     const DfEriEngine::Query queryRS(0, 0, kShellType, lShellType);
//                     const DfEriEngine::CGTO_Pair pairPQ = DfEriEngine::getCGTO_pair(this->m_TlOrbitalInfo,
//                                                                                     iShell, jShell);
//                     const DfEriEngine::CGTO_Pair pairRS = DfEriEngine::getCGTO_pair(this->m_TlOrbitalInfo,
//                                                                                     kShell, lShell);
//                     dfEriEngine.calc(queryPQ, queryRS, pairPQ, pairRS);

// #pragma omp critical (storeK_ID)
//                     {
//                         this->storeKmatrixByIntegralDriven(iShell, iSteps, jShell, jSteps,
//                                                            kShell, kSteps, lShell, lSteps, P, dfEriEngine.WORK, pK);
//                     }
//                 }

                // for debug
                if (this->isOutputTEI_ == true) {
                    std::ofstream fo;
                    fo.open(this->outputTEI_path_.c_str(), std::ios::out | std::ios::app);

                    int index = 0;
                    for (int i = 0; i < iSteps; ++i) {
                        const int I = iShell + i;
                        for (int j = 0; j < jSteps; ++j) {
                            const int J = jShell + j;
                            for (int k = 0; k < kSteps; ++k) {
                                const int K = kShell + k;
                                for (int l = 0; l < lSteps; ++l) {
                                    const int L = lShell + l;
                                    double value = 0.0;
//                                     if (this->isUseNewEngine_ != true) {
                                        value = dfTEI.ERI[index];
//                                     } else {
//                                         value = dfEriEngine.WORK[index];
//                                     }

                                    const std::string output = TlUtils::format("ERI[%d%d%d%d](%5d %5d | %5d %5d) % 16.10f\n",
                                                                               iShellType, jShellType, kShellType, lShellType,
                                                                               I, J, K, L, value);
                                    fo << output;
                                    ++index;
                                }
                            }
                        }
                    }

                    fo.close();
                }

                // HGP
//                 if (this->isOutputHGP_ == true) {
//                     DfHGP::Query eri_IJ;
//                     DfHGP::Query eri_KL;
//                     DfHGP::CGTO_Pair IJ;
//                     DfHGP::CGTO_Pair KL;

//                     eri_IJ.shellType_i = this->m_TlOrbitalInfo.getShellType(iShell);
//                     eri_IJ.shellType_j = this->m_TlOrbitalInfo.getShellType(jShell);
//                     eri_IJ.shellType_i_bar = 0;
//                     eri_IJ.shellType_j_bar = 0;
//                     eri_KL.shellType_i = this->m_TlOrbitalInfo.getShellType(kShell);
//                     eri_KL.shellType_j = this->m_TlOrbitalInfo.getShellType(lShell);
//                     eri_KL.shellType_i_bar = 0;
//                     eri_KL.shellType_j_bar = 0;

//                     const TlPosition A = this->m_TlOrbitalInfo.getPosition(iShell);
//                     const TlPosition B = this->m_TlOrbitalInfo.getPosition(jShell);
//                     IJ.AB = B - A;
//                     IJ.PS.clear();
//                     const int numOfPGTO_I = this->m_TlOrbitalInfo.getCgtoContraction(iShell);
//                     const int numOfPGTO_J = this->m_TlOrbitalInfo.getCgtoContraction(jShell);
//                     for (int i = 0; i < numOfPGTO_I; ++i) {
//                         const int a = this->m_TlOrbitalInfo.getShellType(iShell);
//                         const double coef_i = this->m_TlOrbitalInfo.getCoefficient(iShell, i);
//                         const double exp_i = this->m_TlOrbitalInfo.getExponent(iShell, i);
//                         for (int j = 0; j < numOfPGTO_J; ++j) {
//                             const int b = this->m_TlOrbitalInfo.getShellType(jShell);
//                             const double coef_j = this->m_TlOrbitalInfo.getCoefficient(jShell, j);
//                             const double exp_j = this->m_TlOrbitalInfo.getExponent(jShell, j);

//                             DfHGP::PGTO_Pair pair(a, coef_i, exp_i,
//                                                   b, coef_j, exp_j,
//                                                   A, B);
//                             IJ.PS.push_back(pair);
//                         }
//                     }

//                     const TlPosition C = this->m_TlOrbitalInfo.getPosition(kShell);
//                     const TlPosition D = this->m_TlOrbitalInfo.getPosition(lShell);
//                     KL.AB = D - C;
//                     KL.PS.clear();
//                     const int numOfPGTO_K = this->m_TlOrbitalInfo.getCgtoContraction(kShell);
//                     const int numOfPGTO_L = this->m_TlOrbitalInfo.getCgtoContraction(lShell);
//                     for (int k = 0; k < numOfPGTO_K; ++k) {
//                         const int c = this->m_TlOrbitalInfo.getShellType(kShell);
//                         const double coef_k = this->m_TlOrbitalInfo.getCoefficient(kShell, k);
//                         const double exp_k = this->m_TlOrbitalInfo.getExponent(kShell, k);
//                         for (int l = 0; l < numOfPGTO_L; ++l) {
//                             const int d = this->m_TlOrbitalInfo.getShellType(lShell);
//                             const double coef_l = this->m_TlOrbitalInfo.getCoefficient(lShell, l);
//                             const double exp_l = this->m_TlOrbitalInfo.getExponent(lShell, l);

//                             DfHGP::PGTO_Pair pair(c, coef_k, exp_k,
//                                                   d, coef_l, exp_l,
//                                                   C, D);
//                             KL.PS.push_back(pair);
//                         }
//                     }
                    
//                     std::cerr << TlUtils::format(">>>> HGP (%d %d | %d %d)",
//                                                  iShell, jShell, kShell, lShell)
//                               << std::endl;
//                     dfHGP.calc(eri_IJ, eri_KL, IJ, KL);

//                     {
//                         std::ofstream fo;
//                         fo.open("TEI_HGP.txt", std::ios::out | std::ios::app);
                        
//                         int index = 0;
//                         for (int i = 0; i < iSteps; ++i) {
//                             const int I = iShell + i;
//                             for (int j = 0; j < jSteps; ++j) {
//                                 const int J = jShell + j;
//                                 for (int k = 0; k < kSteps; ++k) {
//                                     const int K = kShell + k;
//                                     for (int l = 0; l < lSteps; ++l) {
//                                         const int L = lShell + l;
//                                         const double value = dfHGP.WORK[index];
                                        
//                                         const std::string output = TlUtils::format("ERI[%d%d%d%d](%5d %5d | %5d %5d) % 16.10f\n",
//                                                                                    iShellType, jShellType, kShellType, lShellType,
//                                                                                    I, J, K, L, value);
//                                         fo << output;
//                                         ++index;
//                                     }
//                                 }
//                             }
//                         }
                        
//                         fo.close();
//                     }
                    
//                 }
            }
        }
    }

}


void DfTwoElectronIntegral::storeKmatrixByIntegralDriven(const int ishell, const int istep, const int jshell, const int jstep,
                                                         const int kshell, const int kstep, const int lshell, const int lstep,
                                                         const TlMatrix& P, const double* pERI, TlMatrix* pKmat)
{
    assert(pKmat != NULL);

    int index = 0;
    for (int i = 0; i < istep; ++i) {
        const int I = ishell + i;

        for (int j = 0; j < jstep; ++j) {
            const int J = jshell + j;
            if (I < J) {
                // bypass
                index += (kstep * lstep);
                continue;
            }

            for (int k = 0; k < kstep; ++k) {
                const int K = kshell + k;
//      if (I < K){
//        // bypass
//        index += lstep;
//        continue;
//      }

                if (jshell == lshell) {
                    if (I < K) {
                        index += lstep;
                        continue;
                    }
                }

                for (int l = 0; l < lstep; ++l) {
                    const int L = lshell + l;
                    if (L > ((I == K) ? J : K)) {
                        // bypass
                        ++index;
                        continue;
                    }

                    const double value = pERI[index];

                    // Eq.1 : (I, K) <= (J, L)
                    if ((I == K) && (J != L)) {
                        (*pKmat)(I, I) -= 2.0*P(J, L) * value;
                    } else {
                        (*pKmat)(I, K) -= P(J, L) * value;
                    }

                    // Eq.2 : (I, L) <= (J, K)
                    if (K != L) {
                        // Eq.1 != Eq.2
                        if ((I == L) && (J != K)) {
                            (*pKmat)(I, I) -= 2.0*P(J, K) * value;
                        } else {
                            (*pKmat)(I, L) -= P(J, K) * value;
                        }
                    }

                    if (I != J) {
                        // (Eq.1, Eq.2) != (Eq.3, Eq.4)

                        // Eq.3 : (J, K) <= (I, L)
                        if ((I != K) || (J != L)) { // Eq.2 != Eq.3 : !((I, L) == (J, K))
                            if ((J == K) && (I != L)) {
                                (*pKmat)(J, J) -= 2.0*P(I, L) * value;
                            } else {
                                (*pKmat)(J, K) -= P(I, L) * value;
                            }
                        }

                        // Eq.4 : (J, L) <= (I, K)
                        if (K != L) {
                            // Eq.3 != Eq.4
                            if ((J == L) && (I != K)) {
                                (*pKmat)(J, J) -= 2.0*P(I, K) * value;
                            } else {
                                (*pKmat)(J, L) -= P(I, K) * value;
                            }
                        }
                    }

                    ++index;
                }
            }

        }
    }
}


// ---------------------------------------------------------------------
// RT method
// ---------------------------------------------------------------------
void DfTwoElectronIntegral::getContractKMatrixByRTmethodCore(const TlSymmetricMatrix& P,
                                                             TlSymmetricMatrix* pK)
{
    assert(pK != NULL);

    // Matrix
    const int nMaxIShell = this->m_TlOrbitalInfo.getNumOfOrbitals();
    pK->resize(nMaxIShell);

    // core
    for (int iShellType = MAX_SHELL_TYPE -1; iShellType >= 0; --iShellType) {

        for (int kShellType = MAX_SHELL_TYPE -1; kShellType >= 0; --kShellType) {
            const int nIKShellPairIndex = 4*iShellType +kShellType;
            const std::vector<IKShellPair> aIKShellPairList = this->getLocalIKShellPairList(nIKShellPairIndex);

            for (int jShellType = MAX_SHELL_TYPE -1; jShellType >= 0; --jShellType) {
                for (int lShellType = MAX_SHELL_TYPE -1; lShellType >= 0; --lShellType) {

                    const unsigned int eriType = DfTEI::getEriType(iShellType, jShellType, kShellType, lShellType);

                    std::vector<IKShellPair>::const_iterator ikShellPairEnd = aIKShellPairList.end();
                    for (std::vector<IKShellPair>::const_iterator ikShellPair = aIKShellPairList.begin();
                            ikShellPair != ikShellPairEnd; ++ikShellPair) {
                        const int iShell = ikShellPair->nIShell;
                        const int kShell = ikShellPair->nKShell;

                        this->RT_Core(eriType, iShell, kShell, P, P, pK);
                    }

                }
            }
        }
    }

    // 後処理
    //(*pK) *= -1.0;

    this->finalize(*pK);
}


void DfTwoElectronIntegral::finalize(TlSymmetricMatrix& K)
{
    // do nothing
}


void DfTwoElectronIntegral::makeShellList()
{
    const int nMaxShell = this->m_TlOrbitalInfo.getNumOfOrbitals();

    int nShell = 0;
    while (nShell < nMaxShell) {
        // nShellType: 0=s, 1=p, 2=d
        const int nShellType = this->m_TlOrbitalInfo.getShellType(nShell);
        const int nSteps = 2 * nShellType +1;

        this->m_ShellList[nShellType].push_back(nShell);

        nShell += nSteps;
    }
}


// this->m_ShellList
void DfTwoElectronIntegral::screenDistributions()
{
    const double dCutoffThreshold = this->m_dCutoffThreshold;

    DfTEI dfTEI(m_dPrimitiveIntegralsCutoffThreshold);

    /// (ij|ij)^(1/2)のi-shell(key)に対する最大値
    std::map<int, double> myuMax;

    DfTEI::ShellPair shellPair;

    // for cutoff report
    std::vector<int> totalNumOfIJ(MAX_SHELL_TYPE * MAX_SHELL_TYPE, 0); // 全IJ-shellの数
    std::vector<int> survivalNumOfIJ(MAX_SHELL_TYPE * MAX_SHELL_TYPE, 0); // 生き残ったIJ-shellの数

    // i-shellのshell-typeによるループ
    for (int iShellType = 0; iShellType < MAX_SHELL_TYPE; ++iShellType) {

        // j-shellのshell-typeによるループ
        for (int jShellType = 0; jShellType < MAX_SHELL_TYPE; ++jShellType) {

            // for (ij|ij) calc
            const unsigned int eriType = DfTEI::getEriType(iShellType, jShellType, iShellType, jShellType);
            int nMaxEriCheck = 0;
            {
                int ijSize = (iShellType*2 +1) * (jShellType*2 +1);
                nMaxEriCheck = ijSize * (ijSize +1) / 2;
            }

            // for cutoff report
            const int ijType = (iShellType * MAX_SHELL_TYPE) + jShellType;

            // i-shellの軌道のループ
            std::vector<int>::const_iterator iShellEnd = this->m_ShellList[iShellType].end();
            for (std::vector<int>::const_iterator iShell = this->m_ShellList[iShellType].begin(); iShell != iShellEnd; ++iShell) {
                // setting variables for i-shell
                const int nIShell = *iShell;
                const TlPosition A = this->m_TlOrbitalInfo.getPosition(nIShell);
                const int nContractionsI = this->m_TlOrbitalInfo.getCgtoContraction(nIShell);
                shellPair.A = A;
                if (myuMax.find(nIShell) == myuMax.end()) {
                    myuMax[nIShell] = 0.0;
                }

                // j-shellの軌道のループ
                std::vector<int>::const_iterator jShellEnd = std::upper_bound(this->m_ShellList[jShellType].begin(),
                                                                              this->m_ShellList[jShellType].end(), nIShell);
                for (std::vector<int>::const_iterator jShell = this->m_ShellList[jShellType].begin(); jShell != jShellEnd; ++jShell) {
                    // setting variables for j-shell
                    const int nJShell = *jShell;
                    // 対称性  assert(nIShell >= nJShell);

                    const TlPosition B = this->m_TlOrbitalInfo.getPosition(nJShell);
                    const int nContractionsJ = this->m_TlOrbitalInfo.getCgtoContraction(nJShell);
                    shellPair.AB = A - B;
                    if (myuMax.find(nJShell) == myuMax.end()) {
                        myuMax[nJShell] = 0.0;
                    }

                    // for primitives
                    const double dRab = A.squareDistanceFrom(B);

                    shellPair.pList.clear();
                    // i-shell primitive loop
                    for (int i = 0; i < nContractionsI; ++i) {
                        const double zeta_a = this->m_TlOrbitalInfo.getExponent(nIShell, i); // i番目のpGTOの軌道指数
                        const double coef_a = this->m_TlOrbitalInfo.getCoefficient(nIShell, i);

                        // j-shell primitive loop
                        for (int j = 0; j < nContractionsJ; ++j) {
                            const double zeta_b = this->m_TlOrbitalInfo.getExponent(nJShell, j);
                            const double coef_b = this->m_TlOrbitalInfo.getCoefficient(nJShell, j);

                            const double zeta = zeta_a + zeta_b;
                            const double inv_zeta = 1.0 / zeta;
                            const double xi = (zeta_a * zeta_b) * inv_zeta;

                            const double xiRab = xi * dRab;
                            const double dKab = CK * std::exp(-xiRab) * inv_zeta * coef_a * coef_b;

                            // cut-off
                            if (std::fabs(dKab) < dCutoffThreshold) {
                                // zetaがTlOrbitalInfoクラスで小さい順にソートされており、
                                // 以降はこれよりも小さいのでj-shell primitive loopを抜ける(break)
                                break;
                            }

                            const TlPosition P = inv_zeta * (zeta_a * A + zeta_b * B);

                            // store ijTable
                            DfTEI::PrimitiveShellPair psp(zeta, dKab, P);
                            shellPair.pList.push_back(psp);
                        }
                    }

                    // for cutoff report
                    ++(totalNumOfIJ[ijType]);

                    if (shellPair.pList.empty() != true) {
                        // primitive listが空で無いものは計算する必要がある(survive)

                        // for cutoff report
                        ++(survivalNumOfIJ[ijType]);

                        // sort
                        // dKabの絶対値の大きい順にソート
                        std::sort(shellPair.pList.begin(), shellPair.pList.end(),
                                  DfTEI::sort_primitiveShellPair_cmp());

                        // keep max |(ij|ij)| value
                        //(this->*pContractEriFunc)(shellPair, shellPair);
                        dfTEI.calc(eriType, shellPair, shellPair);
                        double dMaxValue = 0.0;
                        for (int i = 0; i < nMaxEriCheck; ++i) {
                            dMaxValue = std::max(dMaxValue, std::fabs(dfTEI.ERI[i]));
                        }
                        dMaxValue = std::sqrt(dMaxValue);
                        shellPair.dSchwartzValue = dMaxValue;
                        //  for dMyuMax
                        myuMax[nIShell] = std::max(myuMax[nIShell], dMaxValue);
                        myuMax[nJShell] = std::max(myuMax[nJShell], dMaxValue);

                        // store to shell pair
                        const ShellPairIndex spIndex(nIShell, nJShell);
                        this->m_storedShellPair[spIndex] = shellPair;

                        // store to distribution list
                        this->m_DistributionShellList[nIShell][jShellType].push_back(nJShell);
                        if (nIShell != nJShell) {
                            this->m_DistributionShellList[nJShell][iShellType].push_back(nIShell);
                        }
                    }

                } // end of j-shell loop
            }

        } // end of i-shell loop
    }

    // sort distribution list
    {
        std::map<int, std::map<int, std::vector<int> > >::iterator iShellEnd = this->m_DistributionShellList.end();
        for (std::map<int, std::map<int, std::vector<int> > >::iterator iShell = this->m_DistributionShellList.begin();
                iShell != iShellEnd; ++iShell) {
            for (int nShellType = 0; nShellType < MAX_SHELL_TYPE; ++nShellType) {
                std::sort(iShell->second[nShellType].begin(), iShell->second[nShellType].end(), std::less<int>());
            }
        }
    }

    // cutoff report
    this->logger(" screen distributions cutoff report:\n");
    this->cutoffReport("ss", totalNumOfIJ[0] - survivalNumOfIJ[0], totalNumOfIJ[0]);
    this->cutoffReport("sp", totalNumOfIJ[1] - survivalNumOfIJ[1], totalNumOfIJ[1]);
    this->cutoffReport("sd", totalNumOfIJ[2] - survivalNumOfIJ[2], totalNumOfIJ[2]);
    this->cutoffReport("ps", totalNumOfIJ[3] - survivalNumOfIJ[3], totalNumOfIJ[3]);
    this->cutoffReport("pp", totalNumOfIJ[4] - survivalNumOfIJ[4], totalNumOfIJ[4]);
    this->cutoffReport("pd", totalNumOfIJ[5] - survivalNumOfIJ[5], totalNumOfIJ[5]);
    this->cutoffReport("ds", totalNumOfIJ[6] - survivalNumOfIJ[6], totalNumOfIJ[6]);
    this->cutoffReport("dp", totalNumOfIJ[7] - survivalNumOfIJ[7], totalNumOfIJ[7]);
    this->cutoffReport("dd", totalNumOfIJ[8] - survivalNumOfIJ[8], totalNumOfIJ[8]);
    this->logger("\n");
}


void DfTwoElectronIntegral::cutoffReport(const std::string& shell, const long cutoffCount, const long totalCount)
{
    double rate = 0.0;
    if (totalCount > 0) {
        rate = (double)cutoffCount / (double)totalCount * 100.0;
    }
    this->logger(TlUtils::format(" %s shell %10d / %10d (%6.2f%%)\n",
                                 shell.c_str(), cutoffCount, totalCount, rate));
}


void DfTwoElectronIntegral::getShellList_density_nocut()
{
    this->m_IKShellPairList.clear();
    this->m_IKShellPairList.resize(4 * MAX_SHELL_TYPE +MAX_SHELL_TYPE);

    // i-shell type loop
    for (int iShellType = 0; iShellType < MAX_SHELL_TYPE; ++iShellType) {

        // k-shell type loop
        for (int kShellType = 0; kShellType < MAX_SHELL_TYPE; ++kShellType) {
            const int nIKShellPairIndex = 4*iShellType +kShellType;

            // i-shell loop
            std::vector<int>::const_iterator iShellEnd = this->m_ShellList[iShellType].end();
            for (std::vector<int>::const_iterator iShell = this->m_ShellList[iShellType].begin(); iShell != iShellEnd; ++iShell) {
                const int nIShell = *iShell;

                // k-shell loop
                const std::vector<int> KShellList = this->getShellList(this->m_ShellList[kShellType], nIShell +1);
                std::vector<int>::const_iterator kShellEnd = KShellList.end();
                for (std::vector<int>::const_iterator kShell = KShellList.begin(); kShell != kShellEnd; ++kShell) {
                    const int nKShell = *kShell;

                    IKShellPair ikShellPair(nIShell, nKShell);
                    this->m_IKShellPairList[nIKShellPairIndex].push_back(ikShellPair);
                }
            }
        }
    }
}


// shell listの中からnMaxShellまでのリストを抽出する
// rTargetShellListはソート済みであること
std::vector<int> DfTwoElectronIntegral::getShellList(const std::vector<int>& rTargetShellList,
                                                     const int nMaxShell)
{
    // nMaxShell までのリストを抜き出す
    // lower_bound()
    // ソートされたシーケンス[first,last)に対し、valueより大きいか等しい最初の要素の位置を返します。
    std::vector<int>::const_iterator pBegin = rTargetShellList.begin();
    std::vector<int>::const_iterator pEnd = rTargetShellList.end();
    std::vector<int>::const_iterator pFind = std::lower_bound(pBegin, pEnd, nMaxShell);
    std::vector<int> aAnswer(pBegin, pFind);

    return aAnswer;
}



// shell listの中からnMaxShellまでのリストを抽出する
std::vector<int> DfTwoElectronIntegral::getShellList_distribute(const std::vector<int>& rShellList,
                                                                const int nTargetShell, const int nTargetShellType)
{

    const int nMaxSize = std::min(rShellList.size(), this->m_DistributionShellList[nTargetShell][nTargetShellType].size());
    std::vector<int>::const_iterator pDistBegin = this->m_DistributionShellList[nTargetShell][nTargetShellType].begin();
    std::vector<int>::const_iterator pDistEnd = this->m_DistributionShellList[nTargetShell][nTargetShellType].end();

    // 2つの集合の積
    std::vector<int> aAnswer(nMaxSize);
    std::vector<int>::iterator pBegin = aAnswer.begin();
    std::vector<int>::iterator pEnd = std::set_intersection(rShellList.begin(), rShellList.end(),
                                                            pDistBegin, pDistEnd, pBegin);
    aAnswer.resize(std::distance(pBegin, pEnd));

    // sort
//   struct shell_list_sort_functor {
//     bool operator()(int index1, int index2){
//       const ShellPairIndex spi1(nTargetShell, index1);
//       const ShellPairIndex spi2(nTargetShell, index2);

//       const ShellPair sp1 = this->m_storedShellPair[spi1];
//       const ShellPair sp2 = this->m_storedShellPair[spi2];

//       return (std::fabs(sp1.pList[0].dKab) > std::fabs(sp2.pList[0].dKab));
//     }
//   };

//   std::sort(aAnswer.begin(), aAnswer.end(), shell_list_sort_functor());

    return aAnswer;
}


std::size_t DfTwoElectronIntegral::getShellList_distribute(const std::vector<int>& rShellList,
                                                           const int nTargetShell, const int nTargetShellType,
                                                           index_type* pOutShellList)
{

    std::vector<int>::const_iterator pDistBegin = this->m_DistributionShellList[nTargetShell][nTargetShellType].begin();
    std::vector<int>::const_iterator pDistEnd = this->m_DistributionShellList[nTargetShell][nTargetShellType].end();

    // 2つの集合の積
    index_type* p = std::set_intersection(rShellList.begin(), rShellList.end(),
                                          pDistBegin, pDistEnd,
                                          pOutShellList);

    // sort
//   struct shell_list_sort_functor {
//     bool operator()(int index1, int index2){
//       const ShellPairIndex spi1(nTargetShell, index1);
//       const ShellPairIndex spi2(nTargetShell, index2);

//       const ShellPair sp1 = this->m_storedShellPair[spi1];
//       const ShellPair sp2 = this->m_storedShellPair[spi2];

//       return (std::fabs(sp1.pList[0].dKab) > std::fabs(sp2.pList[0].dKab));
//     }
//   };

//   std::sort(aAnswer.begin(), aAnswer.end(), shell_list_sort_functor());

    return std::distance(pOutShellList, p);
}


DfTEI::ShellPair DfTwoElectronIntegral::getShellPair(const int nIShell, const int nJShell)
{
    const ShellPairIndex spIndex(nIShell, nJShell);
    DfTEI::ShellPair shellPair = this->m_storedShellPair[spIndex];

    if (nIShell < nJShell) {
        shellPair.swap();
    }

    return shellPair;
}


double DfTwoElectronIntegral::getMaxDensityMatrixElement(const int nShell1, const int nStep1,
                                                         const int nShell2, const int nStep2,
                                                         const TlSymmetricMatrix& P)
{
    double dAnswer = 0.0;
    for (int i = 0; i < nStep1; ++i) {
        const int I = nShell1 +i;
        for (int j = 0; j < nStep2; ++j) {
            const int J = nShell2 +j;

            dAnswer = std::max(dAnswer, std::fabs(P(I, J)));
        }
    }

    return dAnswer;
}


// IK shellの組を返す
// 並列化用(シリアルでは全部返す)
std::vector<DfTwoElectronIntegral::IKShellPair> DfTwoElectronIntegral::getLocalIKShellPairList(const int nIKShellPairIndex)
{
    // for serial version
    return (this->m_IKShellPairList[nIKShellPairIndex]);
}


// see E. Schwegler and M. Challacombe, J. Chem. Phys., 105, 2726 (1996).
// "Thresholding contributions to K_ab"
void DfTwoElectronIntegral::getShellList_density()
{
    const double dCutoffThreshold = this->m_dCutoffThreshold * 0.01;

    TlFmt& FmT = TlFmt::getInstance();
    const double coef[6] = { -0.017450254,
                             0.132520568,
                             -0.047915444,
                             0.792267596,
                             -0.583721015,
                             0.697593555
                           };
    // l = 1.0のときのgamma値
    const double gamma1[6] = { 0.03,
                               0.081722098,
                               0.222616709,
                               0.606423483,
                               1.651939972,
                               4.5
                             };
    std::vector<double> gamma(6);
    {
        const double inv_lprime = 1.0 / this->lengthScaleParameter_;
        const double ll = inv_lprime * inv_lprime;
        for (int i = 0; i < 6; ++i) {
            gamma[i] = gamma1[i] * ll;
        }
    }


    // for cutoff report
    long nTotal = 0;
    long nKilled = 0;

    this->m_IKShellPairList.clear();
    this->m_IKShellPairList.resize(4 * MAX_SHELL_TYPE +MAX_SHELL_TYPE);

    // i-shell type loop
    for (int iShellType = 0; iShellType < MAX_SHELL_TYPE; ++iShellType) {
        //const int nISteps = 2 * iShellType +1;

        // k-shell type loop
        for (int kShellType = 0; kShellType < MAX_SHELL_TYPE; ++kShellType) {
            //const int nKSteps = 2 * kShellType +1;
            const int nIKShellPairIndex = 4*iShellType +kShellType;

            // i-shell loop
            std::vector<int>::const_iterator iShellEnd = this->m_ShellList[iShellType].end();
            for (std::vector<int>::const_iterator iShell = this->m_ShellList[iShellType].begin(); iShell != iShellEnd; ++iShell) {
                const int nIShell = *iShell;
                const TlPosition A = this->m_TlOrbitalInfo.getPosition(nIShell);
                const int contractionsA = this->m_TlOrbitalInfo.getCgtoContraction(nIShell);

                // k-shell loop
                const std::vector<int> KShellList = this->getShellList(this->m_ShellList[kShellType], nIShell +1);
                std::vector<int>::const_iterator kShellEnd = KShellList.end();
                for (std::vector<int>::const_iterator kShell = KShellList.begin(); kShell != kShellEnd; ++kShell) {
                    const int nKShell = *kShell;
                    const TlPosition B = this->m_TlOrbitalInfo.getPosition(nKShell);
                    const int contractionsB = this->m_TlOrbitalInfo.getCgtoContraction(nKShell);
                    const double AB2 = A.squareDistanceFrom(B);

                    double judge = 0.0;

                    for (int primitiveA = 0; primitiveA < contractionsA; ++primitiveA) {
                        const double zetaA = this->m_TlOrbitalInfo.getExponent(nIShell, primitiveA);
                        const double coefA = this->m_TlOrbitalInfo.getCoefficient(nIShell, primitiveA);

                        for (int primitiveB = 0; primitiveB < contractionsB; ++primitiveB) {
                            const double zetaB = this->m_TlOrbitalInfo.getExponent(nKShell, primitiveB);
                            const double coefB = this->m_TlOrbitalInfo.getCoefficient(nKShell, primitiveB);

                            const double coefAB = coefA * coefB;
                            const double zetaAB = zetaA * zetaB;
                            const double zetaA_B = zetaA + zetaB;
                            for (int i = 0; i < 6; ++i) {
                                const double zetaAgamma = zetaA * gamma[i];
                                const double zetaBgamma = zetaB * gamma[i];
                                const double param = zetaAB + zetaAgamma + zetaBgamma;
                                const double term1 = CONTRIBUTE_COEF / (std::sqrt(zetaA_B) * param);
                                const double term2 = std::exp(-zetaAB * gamma[i] * AB2 / param);
                                const double T = zetaAB * zetaAB * AB2 / (zetaA_B * param);
                                double term3 = 0.0;
                                FmT.getFmT(0, T, &term3);
                                judge += coefAB * coef[i] * term1 * term2 * term3;
                            }
                        }
                    }

                    if (std::fabs(judge) > dCutoffThreshold) {
                        IKShellPair ikShellPair(nIShell, nKShell);
                        this->m_IKShellPairList[nIKShellPairIndex].push_back(ikShellPair);
                    } else {
                        ++nKilled;
                        //this->logger(TlUtils::format("(%3d|%3d) judge=%e\n", nIShell, nKShell, judge));
                    }

                    ++nTotal;
                }
            }
        }
    }

    // cutoff report
    this->logger(TlUtils::format("density cutoff> %d / %d (% 5.1f%%) cutoff.\n",
                                 nKilled, nTotal, (static_cast<double>(nKilled) / static_cast<double>(nTotal))*100.0));
}


