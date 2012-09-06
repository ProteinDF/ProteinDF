#include <cassert>
#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

#include "DfGridFreeXC.h"
#include "DfFunctional_SVWN.h"
#include "DfFunctional_HFS.h"
#include "TlTime.h"
#include "TlRowVectorMatrix2.h"
#include "TlSystem.h"

const int DfGridFreeXC::MAX_SHELL_TYPE = 2 + 1;

DfGridFreeXC::DfGridFreeXC(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam), pOvpEngines_(NULL),
      orbitalInfo_((*pPdfParam)["coordinates"], (*pPdfParam)["basis_sets"]) {

    this->numOfPQs_ = this->m_nNumOfAOs * (this->m_nNumOfAOs + 1) / 2;

    this->tau_ = 1.0E-10;
    if ((*pPdfParam)["grid_free/CDAM_tau"].getStr().empty() != true) {
        this->tau_ = (*pPdfParam)["grid_free/CDAM_tau"].getDouble();
    }    

    this->epsilon_ = 1.0E-4;
    if ((*pPdfParam)["grid_free/CD_epsilon"].getStr().empty() != true) {
        this->epsilon_ = (*pPdfParam)["grid_free/CD_epsilon"].getDouble();
    }    
}


DfGridFreeXC::~DfGridFreeXC()
{
}


void DfGridFreeXC::buildFxc()
{
    const index_type numOfAOs = this->m_nNumOfAOs;
    const TlSymmetricMatrix P = 0.5 * DfObject::getPpqMatrix<TlSymmetricMatrix>(RUN_RKS, this->m_nIteration -1);
    TlSymmetricMatrix S = DfObject::getSpqMatrix<TlSymmetricMatrix>();

    TlSymmetricMatrix M;
    if (this->XC_engine_ == XC_ENGINE_CD) {
        this->log_.info("begin to create M matrix based on CD.");
        this->getM_byCD(&M);
    } else {
        this->log_.info("begin to create M matrix using 4-center overlap.");
        this->getM(P, &M);
    }
    this->log_.info("begin to generate Fxc using grid-free method.");
    // M.save("M.mat");

    // this->getM_exact(P, &M);
    // M.save("M_exact.mat");

    // tV * S * V == I
    TlMatrix V = DfObject::getXMatrix<TlMatrix>();
    TlMatrix tV = V;
    tV.transpose();
    
    TlSymmetricMatrix M_tilda = tV * M * V;
    // M_tilda.save("Mtilda.mat");
    const int numOfOrthgonalAOs = M_tilda.getNumOfRows();
    
    TlMatrix U;
    TlVector lamda;
    M_tilda.diagonal(&lamda, &U);
    // U.save("U.mat");
    // lamda.save("lamda.vct");
    
    TlSymmetricMatrix F_lamda;
    TlSymmetricMatrix E_lamda;
    this->get_F_lamda(lamda, &F_lamda, &E_lamda);
    
    TlMatrix SVU = S * V * U;
    TlMatrix UVS = SVU;
    UVS.transpose();
    
    // TlSymmetricMatrix Fxc = S * V * U * F_lamda * tU * tV * S;
    TlSymmetricMatrix Fxc = SVU * F_lamda * UVS;
    DfObject::saveFxcMatrix(RUN_RKS, this->m_nIteration, Fxc);
    // TlSymmetricMatrix PE = 2.0 * P * Fxc; 

    TlSymmetricMatrix Exc = SVU * E_lamda * UVS;
    DfObject::saveExcMatrix(RUN_RKS, this->m_nIteration, Exc);
    // this->XC_energy_ = E.dot(P).sum() * 2.0; // rks
    // this->log_.info(TlUtils::format("Exc = %f", this->XC_energy_));

    // {
    //     TlMatrix tSVU = S * V * U;
    //     tSVU.transpose();
    //     tSVU.save("tSVU.mat");

    //     TlMatrix tUtVS = tU * tV * S;
    //     tUtVS.save("tUtVS.mat");
    // }
}


void DfGridFreeXC::getM(const TlSymmetricMatrix& P, TlSymmetricMatrix* pM)
{
    assert(pM != NULL);
    pM->resize(this->m_nNumOfAOs);
    // this->check_.resize(this->m_nNumOfAOs);
    // this->check2_.resize(this->m_nNumOfAOs);

    const TlOrbitalInfo orbitalInfo((*(this->pPdfParam_))["coordinates"],
                                    (*(this->pPdfParam_))["basis_sets"]);

    const TlSparseSymmetricMatrix schwarzTable = this->makeSchwarzTable(orbitalInfo);

    this->createEngines();
    DfTaskCtrl* pDfTaskCtrl = this->getDfTaskCtrlObject();
    // pDfTaskCtrl->setCutoffThreshold(this->cutoffThreshold_);
    pDfTaskCtrl->setCutoffThreshold(0.0);
    pDfTaskCtrl->setCutoffEpsilon_density(0.0);  // cannot use this cutoff
    // pDfTaskCtrl->setCutoffEpsilon_distribution(this->cutoffEpsilon_distribution_);
    pDfTaskCtrl->setCutoffEpsilon_distribution(0.0);

    std::vector<DfTaskCtrl::Task4> taskList;
    bool hasTask = pDfTaskCtrl->getQueue4(orbitalInfo,
                                          schwarzTable,
                                          this->grainSize_,
                                          &taskList,
                                          true);
    while (hasTask == true) {
        this->getM_part(orbitalInfo,
                        taskList,
                        P, pM);
        hasTask = pDfTaskCtrl->getQueue4(orbitalInfo,
                                         schwarzTable,
                                         this->grainSize_,
                                         &taskList);
    }
                    
    this->finalize(pM);

    delete pDfTaskCtrl;
    pDfTaskCtrl = NULL;
    this->destroyEngines();
}


void DfGridFreeXC::createEngines()
{
    assert(this->pOvpEngines_ == NULL);
    
#ifdef _OPENMP
    {
        const int numOfThreads = omp_get_max_threads();
        this->log_.info(TlUtils::format("create OpenMP ERI engine: %d", numOfThreads));
        this->pOvpEngines_ = new DfOverlapEngine[numOfThreads];
    }
#else
    this->pOvpEngines_ = new DfOverlapEngine[1];
#endif // _OPENMP
}


void DfGridFreeXC::destroyEngines()
{
    this->log_.info("delete OpenMP ERI engine");
    if (this->pOvpEngines_ != NULL) {
        delete[] this->pOvpEngines_;
    }
    this->pOvpEngines_ = NULL;
}


DfTaskCtrl* DfGridFreeXC::getDfTaskCtrlObject() const
{
    DfTaskCtrl* pDfTaskCtrl = new DfTaskCtrl(this->pPdfParam_);
    return pDfTaskCtrl;
}


void DfGridFreeXC::finalize(TlSymmetricMatrix* pMtx)
{
    // do nothing
}


TlSparseSymmetricMatrix DfGridFreeXC::makeSchwarzTable(const TlOrbitalInfoObject& orbitalInfo)
{
    this->log_.info("make Schwartz cutoff table: start");
    const index_type maxShellIndex = orbitalInfo.getNumOfOrbitals();
    TlSparseSymmetricMatrix schwarz(maxShellIndex);

    DfOverlapEngine engine;
    // engine.setPrimitiveLevelThreshold(0.0);
    
    for (index_type shellIndexP = 0; shellIndexP < maxShellIndex; ) {
        const int shellTypeP = orbitalInfo.getShellType(shellIndexP);
        const int maxStepsP = 2 * shellTypeP + 1;
        const TlPosition posP = orbitalInfo.getPosition(shellIndexP);
        const DfOverlapEngine::PGTOs pgtosP = DfOverlapEngine::getPGTOs(orbitalInfo, shellIndexP);

        for (index_type shellIndexQ = 0; shellIndexQ < maxShellIndex; ) {
            const int shellTypeQ = orbitalInfo.getShellType(shellIndexQ);
            const int maxStepsQ = 2 * shellTypeQ + 1;
            
            const TlPosition posQ = orbitalInfo.getPosition(shellIndexQ);
            const DfOverlapEngine::PGTOs pgtosQ = DfOverlapEngine::getPGTOs(orbitalInfo, shellIndexQ);
            const DfOverlapEngine::Query query(0, 0, 0, 0,
                                               shellTypeP, shellTypeQ,
                                               shellTypeP, shellTypeQ);

            engine.calc(query, posP, posQ, posP, posQ,
                        pgtosP, pgtosQ, pgtosP, pgtosQ);

            double maxValue = 0.0;
            const int maxIndex = maxStepsP * maxStepsQ;
            for (int index = 0; index < maxIndex; ++index) {
                maxValue = std::max(maxValue, std::fabs(engine.WORK[index]));
            }
            schwarz.set(shellIndexP, shellIndexQ, std::sqrt(maxValue));
            
            shellIndexQ += maxStepsQ;
        }
        shellIndexP += maxStepsP;
    }

    this->log_.info("make Schwartz cutoff table: end");
    return schwarz;
}


void DfGridFreeXC::getM_part(const TlOrbitalInfoObject& orbitalInfo,
                             const std::vector<DfTaskCtrl::Task4>& taskList,
                             const TlMatrixObject& P, TlMatrixObject* pM)
{
    const int taskListSize = taskList.size();
    // const double pairwisePGTO_cutoffThreshold = this->cutoffEpsilon3_;

#pragma omp parallel
    {
        int threadID = 0;
#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif // _OPENMP
        
        // this->pOvpEngines_[threadID].setPrimitiveLevelThreshold(this->cutoffEpsilon3_);

#pragma omp for schedule(runtime)
        for (int i = 0; i < taskListSize; ++i) {
            const index_type shellIndexP = taskList[i].shellIndex1;
            const index_type shellIndexQ = taskList[i].shellIndex2;
            const index_type shellIndexR = taskList[i].shellIndex3;
            const index_type shellIndexS = taskList[i].shellIndex4;
            const int shellTypeP = orbitalInfo.getShellType(shellIndexP);
            const int shellTypeQ = orbitalInfo.getShellType(shellIndexQ);
            const int shellTypeR = orbitalInfo.getShellType(shellIndexR);
            const int shellTypeS = orbitalInfo.getShellType(shellIndexS);
            const int maxStepsP = 2 * shellTypeP + 1;
            const int maxStepsQ = 2 * shellTypeQ + 1;
            const int maxStepsR = 2 * shellTypeR + 1;
            const int maxStepsS = 2 * shellTypeS + 1;
            const TlPosition posP = orbitalInfo.getPosition(shellIndexP);
            const TlPosition posQ = orbitalInfo.getPosition(shellIndexQ);
            const TlPosition posR = orbitalInfo.getPosition(shellIndexR);
            const TlPosition posS = orbitalInfo.getPosition(shellIndexS);
            const DfOverlapEngine::PGTOs pgtosP = DfOverlapEngine::getPGTOs(orbitalInfo, shellIndexP);
            const DfOverlapEngine::PGTOs pgtosQ = DfOverlapEngine::getPGTOs(orbitalInfo, shellIndexQ);
            const DfOverlapEngine::PGTOs pgtosR = DfOverlapEngine::getPGTOs(orbitalInfo, shellIndexR);
            const DfOverlapEngine::PGTOs pgtosS = DfOverlapEngine::getPGTOs(orbitalInfo, shellIndexS);
                        
            const DfOverlapEngine::Query query(0, 0, 0, 0,
                                               shellTypeP, shellTypeQ,
                                               shellTypeR, shellTypeS);

            this->pOvpEngines_[threadID].calc(query,
                                              posP, posQ, posR, posS,
                                              pgtosP, pgtosQ, pgtosR, pgtosS);
                        
            this->storeM(shellIndexP, maxStepsP,
                         shellIndexQ, maxStepsQ,
                         shellIndexR, maxStepsR,
                         shellIndexS, maxStepsS,
                         this->pOvpEngines_[threadID], P, pM);
        }
    }
        
}


void DfGridFreeXC::storeM(const index_type shellIndexP, const int maxStepsP,
                          const index_type shellIndexQ, const int maxStepsQ,
                          const index_type shellIndexR, const int maxStepsR,
                          const index_type shellIndexS, const int maxStepsS,
                          const DfOverlapEngine& engine,
                          const TlMatrixObject& P,
                          TlMatrixObject* pM)
{
    int index = 0;
    for (int i = 0; i < maxStepsP; ++i) {
        const index_type indexP = shellIndexP + i;

        for (int j = 0; j < maxStepsQ; ++j) {
            const index_type indexQ = shellIndexQ + j;
            const double P_pq = P.get(indexP, indexQ);
            
            for (int k = 0; k < maxStepsR; ++k) {
                const index_type indexR = shellIndexR + k;

                for (int l = 0; l < maxStepsS; ++l) {
                    const index_type indexS = shellIndexS + l;
                    const double P_rs = P.get(indexR, indexS);
                    
                    const double value = engine.WORK[index];
                    const index_type maxIndexS = (indexP == indexR) ? indexQ : indexR;
                    
                    if ((indexP >= indexQ) && (maxIndexS >= indexS)) {
                        // Eq.1 : (indexP, indexQ) <= (indexR, indexS)
                        const double coefEq1 = (indexR != indexS) ? 2.0 : 1.0;
                        pM->add(indexP, indexQ, coefEq1 * P_rs * value);
                        
                        // Eq.2 : (indexR, indexS) <= (indexP, indexQ)
                        if ((shellIndexP != shellIndexR) || (shellIndexQ != shellIndexS) || (indexP == indexR)) {
                            if (((indexP + indexQ) != (indexR + indexS)) ||
                                ((indexP * indexQ) != (indexR * indexS))) {
                                // Eq.1の条件と重複しないようにするための措置
                                            
                                const double coefEq2 = (indexP != indexQ) ? 2.0 : 1.0;
                                pM->add(indexR, indexS, coefEq2 * P_pq * value);
                            }
                        }
                    }
                    ++index;
                }
            }
        }
    }
}


void DfGridFreeXC::get_F_lamda(const TlVector lamda,
                               TlSymmetricMatrix* pF_lamda,
                               TlSymmetricMatrix* pE_lamda)
{
    const int dim = lamda.getSize();
    pF_lamda->resize(dim);
    pE_lamda->resize(dim);

    DfFunctional_LDA* pFunc = NULL;
    std::string checkXC = this->m_sXCFunctional;
    if (checkXC == "SVWN") {
        pFunc = new DfFunctional_SVWN();
    } else if (checkXC == "HFS") {
        pFunc = new DfFunctional_HFS();
    } else {
        this->log_.critical(TlUtils::format("not support functional: %s", checkXC.c_str()));
        abort();
    }

    double fv_a = 0.0;
    double fv_b = 0.0;
    for (int i = 0; i < dim; ++i) {
        const double v = lamda.get(i);
        if (v > 1.0E-16) {
            pFunc->getDerivativeFunctional(v, v, &fv_a, &fv_b);
            pF_lamda->set(i, i, fv_a);

            const double f = pFunc->getFunctional(v, v) / (2.0 * v);
            pE_lamda->set(i, i, f);
        }
    }

    delete pFunc;
    pFunc = NULL;
}


void DfGridFreeXC::getM_exact(const TlSymmetricMatrix& P, TlSymmetricMatrix* pM)
{
    assert(pM != NULL);
    TlMatrix M(this->m_nNumOfAOs, this->m_nNumOfAOs);
    pM->resize(this->m_nNumOfAOs);

    DfOverlapEngine engine;
    
    const TlOrbitalInfo orbitalInfo((*(this->pPdfParam_))["coordinates"],
                                    (*(this->pPdfParam_))["basis_sets"]);

    const ShellArrayTable shellArrayTable = this->makeShellArrayTable(orbitalInfo);
    // const ShellPairArrayTable shellPairArrayTable = this->getShellPairArrayTable(shellArrayTable);

    for (int shellTypeP = MAX_SHELL_TYPE -1; shellTypeP >= 0; --shellTypeP) {
        const int maxStepsP = 2 * shellTypeP + 1;
        const ShellArray shellArrayP = shellArrayTable[shellTypeP];
        ShellArray::const_iterator pItEnd = shellArrayP.end();

        for (int shellTypeQ = MAX_SHELL_TYPE -1; shellTypeQ >= 0; --shellTypeQ) {
            const int maxStepsQ = 2 * shellTypeQ + 1;
            const ShellArray shellArrayQ = shellArrayTable[shellTypeQ];
            ShellArray::const_iterator qItEnd = shellArrayQ.end();
            
            for (int shellTypeR = MAX_SHELL_TYPE -1; shellTypeR >= 0; --shellTypeR) {
                const int maxStepsR = 2 * shellTypeR + 1;
                const ShellArray shellArrayR = shellArrayTable[shellTypeR];
                ShellArray::const_iterator rItEnd = shellArrayR.end();

                for (int shellTypeS = MAX_SHELL_TYPE -1; shellTypeS >= 0; --shellTypeS) {
                    const int maxStepsS = 2 * shellTypeS + 1;
                    const ShellArray shellArrayS = shellArrayTable[shellTypeS];
                    ShellArray::const_iterator sItEnd = shellArrayS.end();

                    const DfOverlapEngine::Query query(0, 0, 0, 0,
                                                       shellTypeP, shellTypeQ,
                                                       shellTypeR, shellTypeS);

                    for (ShellArray::const_iterator pIt = shellArrayP.begin(); pIt != pItEnd; ++pIt) {
                        const index_type shellIndexP = *pIt;
                        const TlPosition posP = orbitalInfo.getPosition(shellIndexP);
                        const DfOverlapEngine::PGTOs pgtosP = DfOverlapEngine::getPGTOs(orbitalInfo, shellIndexP);

                        for (ShellArray::const_iterator qIt = shellArrayQ.begin(); qIt != qItEnd; ++qIt) {
                            const index_type shellIndexQ = *qIt;
                            const TlPosition posQ = orbitalInfo.getPosition(shellIndexQ);
                            const DfOverlapEngine::PGTOs pgtosQ = DfOverlapEngine::getPGTOs(orbitalInfo, shellIndexQ);

                            for (ShellArray::const_iterator rIt = shellArrayR.begin(); rIt != rItEnd; ++rIt) {
                                const index_type shellIndexR = *rIt;
                                const TlPosition posR = orbitalInfo.getPosition(shellIndexR);
                                const DfOverlapEngine::PGTOs pgtosR = DfOverlapEngine::getPGTOs(orbitalInfo, shellIndexR);
                                
                                for (ShellArray::const_iterator sIt = shellArrayS.begin(); sIt != sItEnd; ++sIt) {
                                    const index_type shellIndexS = *sIt;
                                    const TlPosition posS = orbitalInfo.getPosition(shellIndexS);
                                    const DfOverlapEngine::PGTOs pgtosS = DfOverlapEngine::getPGTOs(orbitalInfo, shellIndexS);

                                    engine.calc(query,
                                                posP, posQ, posR, posS,
                                                pgtosP, pgtosQ, pgtosR, pgtosS);

                                    int index = 0;
                                    for (int i = 0; i < maxStepsP; ++i) {
                                        const int indexP = shellIndexP + i;

                                        for (int j = 0; j < maxStepsQ; ++j) {
                                            const int indexQ = shellIndexQ + j;
                                            
                                            for (int k = 0; k < maxStepsR; ++k) {
                                                const int indexR = shellIndexR + k;

                                                for (int l = 0; l < maxStepsS; ++l) {
                                                    const int indexS = shellIndexS + l;
                                                    
                                                    const double P_rs = P.get(indexR, indexS);
                                                    const double value = engine.WORK[index];
                                                    M.add(indexP, indexQ, P_rs * value);

                                                    ++index;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    *pM = M;
}


DfGridFreeXC::ShellArrayTable DfGridFreeXC::makeShellArrayTable(const TlOrbitalInfoObject& orbitalInfo)
{
    ShellArrayTable shellArrayTable(MAX_SHELL_TYPE);
    const index_type maxShellIndex = orbitalInfo.getNumOfOrbitals();

    index_type shellIndex = 0;
    while (shellIndex < maxShellIndex) {
        // shellType: 0=s, 1=p, 2=d
        const int shellType = orbitalInfo.getShellType(shellIndex);
        const int steps = 2 * shellType +1;

        shellArrayTable[shellType].push_back(shellIndex);
        
        shellIndex += steps;
    }

    return shellArrayTable;
}


DfGridFreeXC::ShellPairArrayTable DfGridFreeXC::getShellPairArrayTable(const ShellArrayTable& shellArrayTable)
{
    ShellPairArrayTable shellPairArrayTable(MAX_SHELL_TYPE * MAX_SHELL_TYPE);

    for (int shellTypeP = MAX_SHELL_TYPE -1; shellTypeP >= 0; --shellTypeP) {
        const ShellArray& shellArrayP = shellArrayTable[shellTypeP];
        ShellArray::const_iterator pItEnd = shellArrayP.end();

        for (int shellTypeR = MAX_SHELL_TYPE -1; shellTypeR >= 0; --shellTypeR) {
            const ShellArray& shellArrayR = shellArrayTable[shellTypeR];
            ShellArray::const_iterator rItEnd = shellArrayR.end();

            const int shellPairType_PR = shellTypeP * MAX_SHELL_TYPE + shellTypeR;
            for (ShellArray::const_iterator pIt = shellArrayP.begin(); pIt != pItEnd; ++pIt) {
                const index_type indexP = *pIt;

                for (ShellArray::const_iterator rIt = shellArrayR.begin(); rIt != rItEnd; ++rIt) {
                    const index_type indexR = *rIt;

                    if (indexP >= indexR) {
                        ShellPair shellPair(indexP, indexR);
                        shellPairArrayTable[shellPairType_PR].push_back(shellPair);
                    }
                }
            }
        }
    }

    return shellPairArrayTable;
}


void DfGridFreeXC::calcCholeskyVectors_onTheFly()
{
    // timing data
    TlTime CD_all_time;
    TlTime CD_diagonals_time;
    TlTime CD_resizeL_time;
    TlTime CD_pivot_time;
    TlTime CD_ERI_time;
    TlTime CD_calc_time;
    TlTime CD_save_time;

    CD_all_time.start();

    this->createEngines();
    this->initializeCutoffStats();

    // this->log_.info(TlUtils::format("# of PQ dimension: %d", int(this->numOfPQs_)));
    TlSparseSymmetricMatrix schwartzTable(this->m_nNumOfAOs);
    PQ_PairArray I2PQ;
    TlVector d; // 対角成分

    CD_diagonals_time.start();
    this->calcDiagonals(&schwartzTable, &I2PQ, &d);
    CD_diagonals_time.stop();

    this->log_.info(TlUtils::format("# of I~ dimension: %d", int(I2PQ.size())));
    this->saveI2PQ(I2PQ);

    const index_type N = I2PQ.size();
    double error = d.getMaxAbsoluteElement();
    std::vector<TlVector::size_type> pivot(N);
    for (index_type i = 0; i < N; ++i) {
        pivot[i] = i;
    }

    // prepare variables
    const bool isUsingMemManager = this->isEnableMmap_;
    TlRowVectorMatrix2 L(N, 1, 1, 0, isUsingMemManager);
    const double threshold = this->epsilon_;
    this->log_.info(TlUtils::format("Cholesky Decomposition: epsilon=%e", this->epsilon_));

    int progress = 0;
    index_type division =  index_type(N * 0.01);
    L.reserve_cols(division);
    index_type m = 0;
    while (error > threshold) {
#ifdef DEBUG_CD
        this->log_.debug(TlUtils::format("CD progress: %12d/%12d: err=% 8.3e", m, N, error));
#endif //DEBUG_CD

        CD_resizeL_time.start();
        // progress 
        if (m >= progress * division) {
            this->log_.info(TlUtils::format("CD progress: %12d: err=% 8.3e, local mem:%8.1f MB",
                                            m, error, TlSystem::getMaxRSS()));
            ++progress;

            // メモリの確保
            L.reserve_cols(division * progress);
        }
        L.resize(N, m+1);
        CD_resizeL_time.stop();

        // pivot
        CD_pivot_time.start();
        {
            std::vector<TlVector::size_type>::const_iterator it = d.argmax(pivot.begin() + m,
                                                                           pivot.end());
            const index_type i = it - pivot.begin();
            std::swap(pivot[m], pivot[i]);
        }
        CD_pivot_time.stop();

        error = d[pivot[m]];
        const double l_m_pm = std::sqrt(d[pivot[m]]);
        L.set(pivot[m], m, l_m_pm);
        
        const double inv_l_m_pm = 1.0 / l_m_pm;

        // ERI
        CD_ERI_time.start();
        const index_type pivot_m = pivot[m];
        std::vector<double> G_pm;
        const index_type numOf_G_cols = N -(m+1);
        {
            std::vector<index_type> G_col_list(numOf_G_cols);
            for (index_type c = 0; c < numOf_G_cols; ++c) {
                const index_type pivot_i = pivot[m+1 +c]; // from (m+1) to N
                G_col_list[c] = pivot_i;
            }
            G_pm = this->getSuperMatrixElements(pivot_m, G_col_list, I2PQ, schwartzTable);
        }
        assert(G_pm.size() == numOf_G_cols);
        CD_ERI_time.stop();

        // CD calc
        CD_calc_time.start();
        const TlVector L_pm = L.getRowVector(pivot_m);
        std::vector<double> L_xm(numOf_G_cols);
#pragma omp parallel for schedule(runtime)
        for (index_type i = 0; i < numOf_G_cols; ++i) {
            const index_type pivot_i = pivot[m+1 +i]; // from (m+1) to N
            TlVector L_pi = L.getRowVector(pivot_i);
            const double sum_ll = (L_pi.dot(L_pm)).sum();
            const double l_m_pi = (G_pm[i] - sum_ll) * inv_l_m_pm;

#pragma omp atomic
            L_xm[i] += l_m_pi; // for OpenMP
            
#pragma omp atomic
            d[pivot_i] -= l_m_pi * l_m_pi;
        }
        for (index_type i = 0; i < numOf_G_cols; ++i) {
            const index_type pivot_i = pivot[m+1 +i]; // from (m+1) to N
            L.set(pivot_i, m, L_xm[i]);
        }
        CD_calc_time.stop();

        ++m;
    }
    L.resize(N, m);
    this->log_.info(TlUtils::format("Cholesky Vectors: %d", m));

    this->destroyEngines();
    //this->schwartzCutoffReport();

    CD_save_time.start();
    this->saveL(L.getTlMatrix());
    CD_save_time.stop();

    CD_all_time.stop();

    // timing data
    this->log_.info(TlUtils::format("CD all:       %12.2f sec.", CD_all_time.getElapseTime()));
    this->log_.info(TlUtils::format("CD diagonals: %12.2f sec.", CD_diagonals_time.getElapseTime()));
    this->log_.info(TlUtils::format("CD resize L:  %12.2f sec.", CD_resizeL_time.getElapseTime()));
    this->log_.info(TlUtils::format("CD pivot:     %12.2f sec.", CD_pivot_time.getElapseTime()));
    this->log_.info(TlUtils::format("CD ERI:       %12.2f sec.", CD_ERI_time.getElapseTime()));
    this->log_.info(TlUtils::format("CD calc:      %12.2f sec.", CD_calc_time.getElapseTime()));
    this->log_.info(TlUtils::format("CD save:      %12.2f sec.", CD_save_time.getElapseTime()));
}


void DfGridFreeXC::calcDiagonals(TlSparseSymmetricMatrix *pSchwartzTable,
                                 PQ_PairArray *pI2PQ,
                                 TlVector *pDiagonals)
{
    const index_type numOfAOs = this->m_nNumOfAOs;
    assert(numOfAOs == this->orbitalInfo_.getNumOfOrbitals());

    const double tau = this->tau_;
    this->log_.info(TlUtils::format("CDAM tau: %e", tau));
    // this->log_.info(TlUtils::format("primitive GTO quartet threshold: %e", this->cutoffEpsilon3_));

    // initialize
    pI2PQ->clear();
    pI2PQ->reserve(this->numOfPQs_);
    pSchwartzTable->clear();
    pSchwartzTable->resize(numOfAOs);
    TlSparseSymmetricMatrix diagonalMat(numOfAOs);

    // task
    this->log_.info("diagonal calculation: start");
    DfTaskCtrl* pDfTaskCtrl = this->getDfTaskCtrlObject();
    std::vector<DfTaskCtrl::Task2> taskList;
    bool hasTask = pDfTaskCtrl->getQueue2(this->orbitalInfo_,
                                          true,
                                          this->grainSize_,
                                          &taskList, true);
    while (hasTask == true) {
        this->calcDiagonals_kernel(taskList,
                                   pSchwartzTable,
                                   &diagonalMat, pI2PQ);
        hasTask = pDfTaskCtrl->getQueue2(this->orbitalInfo_,
                                         true,
                                         this->grainSize_,
                                         &taskList);
    }
    delete pDfTaskCtrl;
    pDfTaskCtrl = NULL;
    
    // finalize
    this->log_.info("diagonal calculation: finalize");
    // this->finalize_I2PQ(pI2PQ);
    // this->finalize(&diagonalMat);
    // this->finalize(pSchwartzTable);
    std::sort(pI2PQ->begin(), pI2PQ->end());

    // set diagonals
    const index_type numOfI = pI2PQ->size();
    pDiagonals->resize(numOfI);
    for (index_type i = 0; i < numOfI; ++i) {
        const index_type row = (*pI2PQ)[i].index1();
        const index_type col = (*pI2PQ)[i].index2();
        const double value = diagonalMat.get(row, col);
        (*pDiagonals)[i] = value;
    }
}


void DfGridFreeXC::calcDiagonals_kernel(const std::vector<DfTaskCtrl::Task2>& taskList,
                                        TlSparseSymmetricMatrix *pSchwartzTable,
                                        TlSparseSymmetricMatrix *pDiagonalMat,
                                        PQ_PairArray *pI2PQ)
{
    const index_type numOfAOs = this->orbitalInfo_.getNumOfOrbitals();
    pDiagonalMat->resize(numOfAOs);

    const double tau = this->tau_;
    const int taskListSize = taskList.size();
    // const double pairwisePGTO_cutoffThreshold = this->cutoffEpsilon3_;

#pragma omp parallel
    {
        PQ_PairArray local_I2PQ;
        TlSparseSymmetricMatrix local_diagMat(numOfAOs);
        TlSparseSymmetricMatrix local_schwartzTable(numOfAOs);
        int threadID = 0;
#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif // _OPENMP
        // this->pEriEngines_[threadID].setPrimitiveLevelThreshold(this->cutoffEpsilon3_);
        
#pragma omp for schedule(runtime)
        for (int i = 0; i < taskListSize; ++i) {
            const index_type shellIndexP = taskList[i].shellIndex1;
            const index_type shellIndexQ = taskList[i].shellIndex2;
            const int shellTypeP = this->orbitalInfo_.getShellType(shellIndexP);
            const int shellTypeQ = this->orbitalInfo_.getShellType(shellIndexQ);
            const int maxStepsP = 2 * shellTypeP + 1;
            const int maxStepsQ = 2 * shellTypeQ + 1;
            const TlPosition posP = this->orbitalInfo_.getPosition(shellIndexP);
            const TlPosition posQ = this->orbitalInfo_.getPosition(shellIndexQ);
            const DfOverlapEngine::PGTOs pgtosP = DfOverlapEngine::getPGTOs(this->orbitalInfo_, shellIndexP);
            const DfOverlapEngine::PGTOs pgtosQ = DfOverlapEngine::getPGTOs(this->orbitalInfo_, shellIndexQ);

            const DfOverlapEngine::Query query(0, 0, 0, 0,
                                               shellTypeP, shellTypeQ,
                                               shellTypeP, shellTypeQ);
            this->pOvpEngines_[threadID].calc(query,
                                              posP, posQ, posP, posQ,
                                              pgtosP, pgtosQ, pgtosP, pgtosQ);
                
            const int maxStepsPQ = maxStepsP * maxStepsQ;
            double maxValue = 0.0;
            for (int p = 0; p < maxStepsP; ++p) {
                const index_type indexP = shellIndexP + p;
                for (int q = 0; q < maxStepsQ; ++q) {
                    const index_type indexQ = shellIndexQ + q;
                    
                    if ((shellIndexP != shellIndexQ) || (indexP >= indexQ)) {
                        const int pq_index = p * maxStepsQ + q;
                        const int pqpq_index = pq_index * maxStepsPQ + pq_index;
                        
                        const double value = this->pOvpEngines_[threadID].WORK[pqpq_index];

                        // for schwartz
                        maxValue = std::max(maxValue, std::fabs(value));
                        
                        // for I~ to pq table
                        if (std::fabs(value) > tau) {
                            if (value > 0) {
                                local_diagMat.set(indexP, indexQ, value);
                                local_I2PQ.push_back(IndexPair2(indexP, indexQ));
                            } else {
                                this->log_.warn(TlUtils::format("pqpq_value: (%d %d)=% e is spoiled.",
                                                                indexP, indexQ, value));
                            }
                        }
                    }
                }
            }
            local_schwartzTable.set(shellIndexP, shellIndexQ, std::sqrt(maxValue));
        }

        // add up
#ifdef _OPENMP
        {
            const int numOfThreads = omp_get_num_threads();
            for (int i = 0; i < numOfThreads; ++i) {
                if (threadID == i) {
                    pI2PQ->insert(pI2PQ->end(),
                                  local_I2PQ.begin(), local_I2PQ.end());
                    pDiagonalMat->merge(local_diagMat);
                    pSchwartzTable->merge(local_schwartzTable);
                }
#pragma omp barrier                
            }
        }
#else
        {
            *pI2PQ = local_I2PQ;
            *pDiagonalMat = local_diagMat;
            *pSchwartzTable = local_schwartzTable;
        }
#endif // _OPENMP
    }
}


void DfGridFreeXC::saveI2PQ(const PQ_PairArray& I2PQ) 
{
    std::string filepath = this->getI2pqVtrPath();
    std::ofstream ofs;
    ofs.open(filepath.c_str(), std::ofstream::out | std::ofstream::binary);

    const std::size_t size = I2PQ.size();
    ofs.write(reinterpret_cast<const char*>(&size), sizeof(std::size_t));
    for (std::size_t i = 0; i < size; ++i) {
        const index_type index1 = I2PQ[i].index1();
        const index_type index2 = I2PQ[i].index2();
        ofs.write(reinterpret_cast<const char*>(&index1), sizeof(index_type));
        ofs.write(reinterpret_cast<const char*>(&index2), sizeof(index_type));
    }

    ofs.close();
}


void DfGridFreeXC::saveL(const TlMatrix& L)
{
    L.save("GF_L.mat");
}


std::vector<double>
DfGridFreeXC::getSuperMatrixElements(const index_type G_row,
                                     const std::vector<index_type>& G_col_list,
                                     const PQ_PairArray& I2PQ,
                                     const TlSparseSymmetricMatrix& schwartzTable)
{
    this->elements_cache_.clear();

    const std::vector<IndexPair4> calcList = this->getCalcList(G_row, G_col_list, I2PQ);
    this->calcElements(calcList, schwartzTable);
    const std::vector<double> answer = this->setElements(G_row, G_col_list, I2PQ);

    return answer;
}


std::vector<DfGridFreeXC::IndexPair4> 
DfGridFreeXC::getCalcList(const index_type G_row,
                          const std::vector<index_type>& G_col_list,
                          const PQ_PairArray& I2PQ)
{
    std::set<IndexPair4> calcSet;

    const index_type indexP = I2PQ[G_row].index1();
    const index_type indexQ = I2PQ[G_row].index2();
    const index_type shellIndexP = this->orbitalInfo_.getShellIndex(indexP);
    const index_type shellIndexQ = this->orbitalInfo_.getShellIndex(indexQ);

    const index_type numOf_G_cols = G_col_list.size();
#pragma omp parallel for schedule(runtime)
    for (index_type i = 0; i < numOf_G_cols; ++i) {
        const index_type G_col = G_col_list[i];

        const index_type indexR = I2PQ[G_col].index1();
        const index_type indexS = I2PQ[G_col].index2();
        const index_type shellIndexR = this->orbitalInfo_.getShellIndex(indexR);
        const index_type shellIndexS = this->orbitalInfo_.getShellIndex(indexS);

        IndexPair4 indexPair4(shellIndexP, shellIndexQ, shellIndexR, shellIndexS);
#pragma omp critical(DfGridFreeXC__getCalcList) 
        {
            calcSet.insert(indexPair4);
        }
    }

    std::vector<IndexPair4> calcList(calcSet.size());
    std::copy(calcSet.begin(), calcSet.end(), calcList.begin());

    return calcList;
}


void DfGridFreeXC::calcElements(const std::vector<IndexPair4>& calcList,
                                const TlSparseSymmetricMatrix& schwartzTable) 
{
    const int maxShellType = this->orbitalInfo_.getMaxShellType();
    const double threshold = this->tau_;
    // const double pairwisePGTO_cutoffThreshold = this->cutoffEpsilon3_;

    const int numOfList = calcList.size();
#pragma omp parallel
    {
        int threadID = 0;
        ElementsCacheType local_cache;

#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif // _OPENMP
        // this->pOvpEngines_[threadID].setPrimitiveLevelThreshold(this->cutoffEpsilon3_);
        
#pragma omp for schedule(runtime)
        for (int i = 0; i < numOfList; ++i) {
            const index_type shellIndexP = calcList[i].index1();
            const index_type shellIndexQ = calcList[i].index2();
            const index_type shellIndexR = calcList[i].index3();
            const index_type shellIndexS = calcList[i].index4();
            
            const int shellTypeP = this->orbitalInfo_.getShellType(shellIndexP);
            const int shellTypeQ = this->orbitalInfo_.getShellType(shellIndexQ);
            const int shellTypeR = this->orbitalInfo_.getShellType(shellIndexR);
            const int shellTypeS = this->orbitalInfo_.getShellType(shellIndexS);
            
            const int shellQuartetType =
                ((shellTypeP * maxShellType + shellTypeQ) * maxShellType + shellTypeP) * maxShellType + shellTypeQ;
            const bool isAlive = this->isAliveBySchwartzCutoff(shellIndexP, shellIndexQ,
                                                               shellIndexR, shellIndexS,
                                                               shellQuartetType,
                                                               schwartzTable,
                                                               threshold);
            if (isAlive == true) {
                const int maxStepsP = 2 * shellTypeP + 1;
                const int maxStepsQ = 2 * shellTypeQ + 1;
                const int maxStepsR = 2 * shellTypeR + 1;
                const int maxStepsS = 2 * shellTypeS + 1;
                
                const DfOverlapEngine::Query query(0, 0, 0, 0,
                                                   shellTypeP, shellTypeQ, 
                                                   shellTypeR, shellTypeS);
                const TlPosition posP = this->orbitalInfo_.getPosition(shellIndexP);
                const TlPosition posQ = this->orbitalInfo_.getPosition(shellIndexQ);
                const TlPosition posR = this->orbitalInfo_.getPosition(shellIndexR);
                const TlPosition posS = this->orbitalInfo_.getPosition(shellIndexS);
                const DfOverlapEngine::PGTOs pgtosP = DfOverlapEngine::getPGTOs(this->orbitalInfo_, shellIndexP);
                const DfOverlapEngine::PGTOs pgtosQ = DfOverlapEngine::getPGTOs(this->orbitalInfo_, shellIndexQ);
                const DfOverlapEngine::PGTOs pgtosR = DfOverlapEngine::getPGTOs(this->orbitalInfo_, shellIndexR);
                const DfOverlapEngine::PGTOs pgtosS = DfOverlapEngine::getPGTOs(this->orbitalInfo_, shellIndexS);
                
                this->pOvpEngines_[threadID].calc(query,
                                                  posP, posQ, posR, posS,
                                                  pgtosP, pgtosQ, pgtosR, pgtosS);
                
                const int steps = maxStepsP * maxStepsQ * maxStepsR * maxStepsS;
                std::vector<double> buf(steps);
                std::copy(this->pOvpEngines_[threadID].WORK, this->pOvpEngines_[threadID].WORK + steps,
                          buf.begin());

                local_cache[calcList[i]] = buf;
            }
        }

        // merge cache
#pragma omp critical(DfGridFreeXC__calcERIs)
        {
            this->elements_cache_.insert(local_cache.begin(), local_cache.end());
        }
    }
}


bool DfGridFreeXC::isAliveBySchwartzCutoff(const index_type shellIndexP,
                                           const index_type shellIndexQ,
                                           const index_type shellIndexR,
                                           const index_type shellIndexS,
                                           const int shellQuartetType,
                                           const TlSparseSymmetricMatrix& schwarzTable,
                                           const double threshold)
{
    bool answer = false;

    const double sqrt_pqpq = schwarzTable.get(shellIndexP, shellIndexQ);
    const double sqrt_rsrs = schwarzTable.get(shellIndexR, shellIndexS);

    if ((sqrt_pqpq * sqrt_rsrs) >= threshold) {
        answer = true;

#pragma omp atomic
        ++(this->cutoffAlive_schwartz_[shellQuartetType]);
    }

#pragma omp atomic
    ++(this->cutoffAll_schwartz_[shellQuartetType]);

    return answer;
}


void DfGridFreeXC::initializeCutoffStats()
{
    // clear cutoff stats
    const int maxShellType = this->orbitalInfo_.getMaxShellType();
    const int numOfShellPairType = maxShellType* maxShellType;
    const int numOfShellQuartetType = numOfShellPairType * numOfShellPairType;
    this->cutoffAll_schwartz_.clear();
    this->cutoffAlive_schwartz_.clear();
    this->cutoffAll_schwartz_.resize(numOfShellQuartetType, 0);
    this->cutoffAlive_schwartz_.resize(numOfShellQuartetType, 0);
}


void DfGridFreeXC::schwartzCutoffReport()
{
    static const char typeStr4[][5] = {
        "SSSS", "SSSP", "SSSD", "SSPS", "SSPP", "SSPD", "SSDS", "SSDP", "SSDD",
        "SPSS", "SPSP", "SPSD", "SPPS", "SPPP", "SPPD", "SPDS", "SPDP", "SPDD",
        "SDSS", "SDSP", "SDSD", "SDPS", "SDPP", "SDPD", "SDDS", "SDDP", "SDDD",
        "PSSS", "PSSP", "PSSD", "PSPS", "PSPP", "PSPD", "PSDS", "PSDP", "PSDD",
        "PPSS", "PPSP", "PPSD", "PPPS", "PPPP", "PPPD", "PPDS", "PPDP", "PPDD",
        "PDSS", "PDSP", "PDSD", "PDPS", "PDPP", "PDPD", "PDDS", "PDDP", "PDDD",
        "DSSS", "DSSP", "DSSD", "DSPS", "DSPP", "DSPD", "DSDS", "DSDP", "DSDD",
        "DPSS", "DPSP", "DPSD", "DPPS", "DPPP", "DPPD", "DPDS", "DPDP", "DPDD",
        "DDSS", "DDSP", "DDSD", "DDPS", "DDPP", "DDPD", "DDDS", "DDDP", "DDDD",
    };
    const int maxShellType = this->orbitalInfo_.getMaxShellType();

    // cutoff report for schwarz
    bool hasCutoffSchwarz = false;
    for (int shellTypeA = 0; shellTypeA < maxShellType; ++shellTypeA) {
        for (int shellTypeB = 0; shellTypeB < maxShellType; ++shellTypeB) {
            const int shellTypeAB = shellTypeA * maxShellType + shellTypeB;
            for (int shellTypeC = 0; shellTypeC < maxShellType; ++shellTypeC) {
                const int shellTypeABC = shellTypeAB * maxShellType + shellTypeC;
                for (int shellTypeD = 0; shellTypeD < maxShellType; ++shellTypeD) {
                    const int shellTypeABCD = shellTypeABC * maxShellType + shellTypeD;
                    if (this->cutoffAll_schwartz_[shellTypeABCD] != 0) {
                        hasCutoffSchwarz = true;
                        break;
                    }
                }
            }
        }
    }
    if (hasCutoffSchwarz == true) {
        this->log_.info("schwarz cutoff report");
        this->log_.info(TlUtils::format("threshold: %e", this->tau_));
        this->log_.info("type: alive / all (ratio)");
        for (int shellTypeA = 0; shellTypeA < maxShellType; ++shellTypeA) {
            for (int shellTypeB = 0; shellTypeB < maxShellType; ++shellTypeB) {
                const int shellTypeAB = shellTypeA * maxShellType + shellTypeB;
                for (int shellTypeC = 0; shellTypeC < maxShellType; ++shellTypeC) {
                    const int shellTypeABC = shellTypeAB * maxShellType + shellTypeC;
                    for (int shellTypeD = 0; shellTypeD < maxShellType; ++shellTypeD) {
                        const int shellTypeABCD = shellTypeABC * maxShellType + shellTypeD;
                        
                        if (this->cutoffAll_schwartz_[shellTypeABCD] > 0) {
                            const double ratio = (double)this->cutoffAlive_schwartz_[shellTypeABCD]
                                / (double)this->cutoffAll_schwartz_[shellTypeABCD]
                                * 100.0;
                            this->log_.info(TlUtils::format(" %4s: %12ld / %12ld (%6.2f%%)",
                                                            typeStr4[shellTypeABCD],
                                                            this->cutoffAlive_schwartz_[shellTypeABCD],
                                                            this->cutoffAll_schwartz_[shellTypeABCD],
                                                            ratio));
                        }
                    }
                }
            }
        }
    }
}


std::vector<double>
DfGridFreeXC::setElements(const index_type G_row,
                          const std::vector<index_type> G_col_list,
                          const PQ_PairArray& I2PQ)
{
    const index_type indexP_orig = I2PQ[G_row].index1();
    const index_type indexQ_orig = I2PQ[G_row].index2();
    const index_type shellIndexP_orig = this->orbitalInfo_.getShellIndex(indexP_orig);
    const index_type shellIndexQ_orig = this->orbitalInfo_.getShellIndex(indexQ_orig);

    const index_type numOf_G_cols = G_col_list.size();
    std::vector<double> answer(numOf_G_cols);
#pragma omp parallel for
    for (index_type i = 0; i < numOf_G_cols; ++i) {
        index_type indexP = indexP_orig;
        index_type indexQ = indexQ_orig;
        index_type shellIndexP = shellIndexP_orig;
        index_type shellIndexQ = shellIndexQ_orig;

        const index_type G_col = G_col_list[i];

        index_type indexR = I2PQ[G_col].index1();
        index_type indexS = I2PQ[G_col].index2();
        index_type shellIndexR = this->orbitalInfo_.getShellIndex(indexR);
        index_type shellIndexS = this->orbitalInfo_.getShellIndex(indexS);

        // swap
        if (shellIndexP > shellIndexQ) {
            std::swap(shellIndexP, shellIndexQ);
            std::swap(indexP, indexQ);
        }
        if (shellIndexR > shellIndexS) {
            std::swap(shellIndexR, shellIndexS);
            std::swap(indexR, indexS);
        }
        {
            IndexPair2 pq(shellIndexP, shellIndexQ);
            IndexPair2 rs(shellIndexR, shellIndexS);
            if (pq.index() > rs.index()) {
                std::swap(shellIndexP, shellIndexR);
                std::swap(indexP, indexR);
                std::swap(shellIndexQ, shellIndexS);
                std::swap(indexQ, indexS);
            }
        }

        std::vector<double> values;
#pragma omp critical(DfGridFreeXC__setElements)
        {
            values = this->elements_cache_[IndexPair4(shellIndexP, shellIndexQ,
                                                      shellIndexR, shellIndexS)];
        }

        if (values.empty() != true) {
            const int basisTypeP = indexP - shellIndexP;
            const int basisTypeQ = indexQ - shellIndexQ;
            const int basisTypeR = indexR - shellIndexR;
            const int basisTypeS = indexS - shellIndexS;
            
            const int shellTypeP = this->orbitalInfo_.getShellType(shellIndexP);
            const int shellTypeQ = this->orbitalInfo_.getShellType(shellIndexQ);
            const int shellTypeR = this->orbitalInfo_.getShellType(shellIndexR);
            const int shellTypeS = this->orbitalInfo_.getShellType(shellIndexS);
            const int maxStepsP = 2 * shellTypeP + 1;
            const int maxStepsQ = 2 * shellTypeQ + 1;
            const int maxStepsR = 2 * shellTypeR + 1;
            const int maxStepsS = 2 * shellTypeS + 1;
            
            const int index = ((basisTypeP * maxStepsQ + basisTypeQ) * maxStepsR + basisTypeR) * maxStepsS + basisTypeS;
            assert(values.size() > index);
            
#pragma omp atomic
            answer[i] += values.at(index);
        }
    }

    return answer;
}


void DfGridFreeXC::getM_byCD(TlSymmetricMatrix* pM)
{
    pM->resize(this->m_nNumOfAOs);
    const TlSymmetricMatrix P = 0.5 * this->getPMatrix();

    // cholesky vector
    TlMatrix L = this->getL();
    const index_type numOfCBs = L.getNumOfCols();

    const PQ_PairArray I2PQ = this->getI2PQ();
    index_type start_CholeskyBasis = 0;
    index_type end_CholeskyBasis = 0;
    this->divideCholeskyBasis(numOfCBs, &start_CholeskyBasis, &end_CholeskyBasis);
    for (index_type I = start_CholeskyBasis; I < end_CholeskyBasis; ++I) {
        TlSymmetricMatrix LI = this->getCholeskyVector(L.getColVector(I), I2PQ);
        assert(LI.getNumOfRows() == this->m_nNumOfAOs);
        assert(LI.getNumOfCols() == this->m_nNumOfAOs);
        
        TlMatrix QI = LI;
        QI.dot(P);
        const double qi = QI.sum();

        *pM += qi*LI;
    }

    // this->finalize(pM);
}


TlSymmetricMatrix DfGridFreeXC::getPMatrix()
{
    TlSymmetricMatrix P = this->getPpqMatrix<TlSymmetricMatrix>(RUN_RKS, this->m_nIteration -1);
    return P;
}


TlMatrix DfGridFreeXC::getL()
{
    //TlMatrix L = DfObject::getLMatrix<TlMatrix>();
    TlMatrix L;
    L.load("GF_L.mat");

    return L;
}


DfGridFreeXC::PQ_PairArray DfGridFreeXC::getI2PQ()
{
    std::string filepath = this->getI2pqVtrPath();
    std::ifstream ifs;
    ifs.open(filepath.c_str(), std::ofstream::in | std::ofstream::binary);
    if (ifs.fail()) {
        abort();
    }

    std::size_t size = 0;
    ifs.read(reinterpret_cast<char*>(&size), sizeof(std::size_t));

    PQ_PairArray answer(size);
    index_type shellIndex1 = 0;
    index_type shellIndex2 = 0;
    for (std::size_t i = 0; i < size; ++i) {
        ifs.read(reinterpret_cast<char*>(&shellIndex1), sizeof(index_type));
        ifs.read(reinterpret_cast<char*>(&shellIndex2), sizeof(index_type));
        answer[i] = IndexPair2(shellIndex1, shellIndex2);
    }

    ifs.close();
    return answer;
}


void DfGridFreeXC::divideCholeskyBasis(const index_type numOfCBs,
                                       index_type *pStart, index_type *pEnd)
{
    *pStart = 0;
    *pEnd = numOfCBs;
}


TlSymmetricMatrix DfGridFreeXC::getCholeskyVector(const TlVector& L_col,
                                                  const PQ_PairArray& I2PQ)
{
    const index_type numOfItilde = L_col.getSize();
    TlSymmetricMatrix answer(this->m_nNumOfAOs);
    for (index_type i = 0; i < numOfItilde; ++i) {
        answer.set(I2PQ[i].index1(),
                   I2PQ[i].index2(),
                   L_col[i]);
    }

    return answer;
}

