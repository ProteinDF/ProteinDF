#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

#include "DfEriX.h"
#include "TlMath.h"
#include "TlFmt.h"

const int DfEriX::MAX_SHELL_TYPE = 2 + 1; // (s=0, p, d)
const int DfEriX::FORCE_K_BUFFER_SIZE = 3 * 5 * 5 * 5; // (xyz) * 5d * 5d * 5d

DfEriX::DfEriX(TlSerializeData* pPdfParam) 
    : DfObject(pPdfParam), pEriEngines_(NULL)
{
    this->lengthScaleParameter_ = 1.0;
    if ((*pPdfParam)["length_scale_parameter"].getStr() != "") {
        this->lengthScaleParameter_ = (*pPdfParam)["length_scale_parameter"].getDouble();
    }
    
    this->cutoffThreshold_ = 1.0E-10;
    if ((*pPdfParam)["cut-value"].getStr().empty() != true) {
        this->cutoffThreshold_ = (*pPdfParam)["cut-value"].getDouble();
    }    

    this->cutoffEpsilon_density_ = this->cutoffThreshold_;
    if ((*pPdfParam)["cutoff_density"].getStr().empty() != true) {
        this->cutoffEpsilon_density_ = (*pPdfParam)["cutoff_density"].getDouble();
    }    

    this->cutoffEpsilon_distribution_ = this->cutoffThreshold_;
    if ((*pPdfParam)["cutoff_distribution"].getStr().empty() != true) {
        this->cutoffEpsilon_distribution_ = (*pPdfParam)["cutoff_distribution"].getDouble();
    }    
    // this->cutoffThreshold_ = 1.0E-10;
    // if ((*pPdfParam)["cut-value"].getStr().empty() != true) {
    //     this->cutoffThreshold_ = (*pPdfParam)["cut-value"].getDouble();
    // }    

    // this->cutoffEpsilon1_ = this->cutoffThreshold_ * 0.01;
    // if ((*pPdfParam)["cutoff_epsilon1"].getStr().empty() != true) {
    //     this->cutoffEpsilon1_ = (*pPdfParam)["cutoff_epsilon1"].getDouble();
    // }    

    // this->cutoffEpsilon2_ = this->cutoffThreshold_;
    // if ((*pPdfParam)["cutoff_epsilon2"].getStr().empty() != true) {
    //     this->cutoffEpsilon2_ = (*pPdfParam)["cutoff_epsilon2"].getDouble();
    // }    

    this->cutoffEpsilon3_ = this->cutoffThreshold_ * 0.01;
    if ((*pPdfParam)["cutoff_epsilon3"].getStr().empty() != true) {
        this->cutoffEpsilon3_ = (*pPdfParam)["cutoff_epsilon3"].getDouble();
    }    

    // debug ===========================================================
    this->isDebugOutJ_ = false;
    if ((*pPdfParam)["debug_out_J"].getStr().empty() != true) {
        this->isDebugOutJ_ = (*pPdfParam)["debug_out_J"].getBoolean();
    }

    this->isDebugOutK_ = false;
    if ((*pPdfParam)["debug_out_K"].getStr().empty() != true) {
        this->isDebugOutK_ = (*pPdfParam)["debug_out_K"].getBoolean();
    }

    this->isDebugExactJ_ = false;
    if ((*pPdfParam)["debug_exact_J"].getStr().empty() != true) {
        this->isDebugExactJ_ = (*pPdfParam)["debug_exact_J"].getBoolean();
    }
    this->isDebugExactK_ = false;
    if ((*pPdfParam)["debug_exact_K"].getStr().empty() != true) {
        this->isDebugExactK_ = (*pPdfParam)["debug_exact_K"].getBoolean();
    }
}


DfEriX::~DfEriX()
{
}


void DfEriX::createEngines()
{
    assert(this->pEriEngines_ == NULL);
    
#ifdef _OPENMP
    {
        const int numOfThreads = omp_get_max_threads();
        this->log_.info(TlUtils::format("create OpenMP ERI engine: %d", numOfThreads));
        this->pEriEngines_ = new DfEriEngine[numOfThreads];
    }
#else
    this->pEriEngines_ = new DfEriEngine[1];
#endif // _OPENMP
}


double DfEriX::destroyEngines()
{
    this->log_.info("delete OpenMP ERI engine");
    
    // sum up elapse time
    double time_ERI_all = 0.0;
    int numOfThreads = 1;
#ifdef _OPENMP
    {
        numOfThreads = omp_get_max_threads();
    }
#endif // _OPENMP
    for (int thread = 0; thread < numOfThreads; ++thread) {
        time_ERI_all += this->pEriEngines_[thread].getElapseCalcTime();
    }

    if (this->pEriEngines_ != NULL) {
        delete[] this->pEriEngines_;
    }
    this->pEriEngines_ = NULL;

    return time_ERI_all;
}


DfTaskCtrl* DfEriX::getDfTaskCtrlObject() const
{
    DfTaskCtrl* pDfTaskCtrl = new DfTaskCtrl(this->pPdfParam_);
    return pDfTaskCtrl;
}


void DfEriX::finalize(TlMatrix* pMtx)
{
    // do nothing
}


void DfEriX::finalize(TlSymmetricMatrix* pMtx)
{
    // do nothing
}


void DfEriX::finalize(TlVector* pVct)
{
    // do nothing
}


void DfEriX::getJ(const TlSymmetricMatrix& P, TlVector* pRho)
{
    assert(pRho != NULL);
    this->elapsetime_store_ = 0.0;
    TlTime time_all;
    time_all.start();
    
    // カットオフ値の設定
    const double maxDeltaP = P.getMaxAbsoluteElement();
    if (maxDeltaP < 1.0) {
        this->cutoffThreshold_ /= std::fabs(maxDeltaP);
        this->log_.info(TlUtils::format(" new cutoff threshold = % e\n",
                                        this->cutoffThreshold_));
    }

    const TlOrbitalInfo orbitalInfo((*(this->pPdfParam_))["coordinates"],
                                    (*(this->pPdfParam_))["basis_sets"]);
    const TlOrbitalInfo_Density orbitalInfo_Density((*(this->pPdfParam_))["coordinates"],
                                                    (*(this->pPdfParam_))["basis_sets_j"]);
    const ShellArrayTable shellArrayTable_Density = this->makeShellArrayTable(orbitalInfo_Density);

    pRho->resize(this->m_nNumOfAux);
    pRho->zeroClear();

    this->createEngines();
    DfTaskCtrl* pDfTaskCtrl = this->getDfTaskCtrlObject();
    pDfTaskCtrl->setCutoffThreshold(this->cutoffThreshold_);

    std::vector<DfTaskCtrl::Task2> taskList;
    bool hasTask = pDfTaskCtrl->getQueue2(orbitalInfo,
                                          true,
                                          this->grainSize_, &taskList, true);
    while (hasTask == true) {
        this->getJ_part(orbitalInfo,
                        orbitalInfo_Density,
                        shellArrayTable_Density,
                        taskList,
                        P, pRho);

        hasTask = pDfTaskCtrl->getQueue2(orbitalInfo,
                                         true,
                                         this->grainSize_, &taskList);
    }

    this->finalize(pRho);

    delete pDfTaskCtrl;
    pDfTaskCtrl = NULL;
    const double time_ERI_all = this->destroyEngines();

    // statics report
    time_all.stop();
    {
        this->log_.info(TlUtils::format("all time:  %16.1f sec.", time_all.getElapseTime()));
        this->log_.info(TlUtils::format(" store:    %16.1f sec.", this->elapsetime_store_));
        this->log_.info(TlUtils::format(" ERI time: %16.1f sec.", time_ERI_all));
    }
}


void DfEriX::getJ_part(const TlOrbitalInfo& orbitalInfo,
                       const TlOrbitalInfo_Density& orbitalInfo_Density,
                       const ShellArrayTable& shellArrayTable_Density,
                       const std::vector<DfTaskCtrl::Task2>& taskList,
                       const TlMatrixObject& P, TlVector* pRho)
{
    const int taskListSize = taskList.size();
    const double pairwisePGTO_cutoffThreshold = this->cutoffEpsilon3_;
    
#pragma omp parallel
    {
        TlTime time_store;
        int threadID = 0;
#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif // _OPENMP

        this->pEriEngines_[threadID].setPrimitiveLevelThreshold(this->cutoffEpsilon3_);

#pragma omp for schedule(runtime)
        for (int i = 0; i < taskListSize; ++i) {
            const index_type shellIndexP = taskList[i].shellIndex1;
            const index_type shellIndexQ = taskList[i].shellIndex2;
            const int shellTypeP = orbitalInfo.getShellType(shellIndexP);
            const int shellTypeQ = orbitalInfo.getShellType(shellIndexQ);
            
            const int maxStepsP = 2 * shellTypeP + 1;
            const int maxStepsQ = 2 * shellTypeQ + 1;
            const TlPosition posP = orbitalInfo.getPosition(shellIndexP);
            const TlPosition posQ = orbitalInfo.getPosition(shellIndexQ);
            const DfEriEngine::CGTO_Pair PQ = this->pEriEngines_[threadID].getCGTO_pair(orbitalInfo,
                                                                                        shellIndexP,
                                                                                        shellIndexQ,
                                                                                        pairwisePGTO_cutoffThreshold);
            const DfEriEngine::Query queryPQ(0, 0, shellTypeP, shellTypeQ);

            for (int shellTypeR = DfEriX::MAX_SHELL_TYPE -1; shellTypeR >= 0; --shellTypeR) {
                const int maxStepsR = 2 * shellTypeR + 1;
            
                const int shellTypeS = 0;
                const DfEriEngine::Query queryRS(0, 0, shellTypeR, shellTypeS);
                
                const std::size_t numOfShellArrayR = shellArrayTable_Density[shellTypeR].size();
                for (std::size_t indexR = 0; indexR < numOfShellArrayR; ++indexR) {
                    const index_type shellIndexR = shellArrayTable_Density[shellTypeR][indexR];
                
                    const DfEriEngine::CGTO_Pair RS = 
                        this->pEriEngines_[threadID].getCGTO_pair(orbitalInfo_Density,
                                                                  shellIndexR,
                                                                  -1,
                                                                  pairwisePGTO_cutoffThreshold);
                    this->pEriEngines_[threadID].calc(queryPQ, queryRS, PQ, RS);
                
                    int index = 0;
                    for (int i = 0; i < maxStepsP; ++i) {
                        const index_type indexP = shellIndexP + i;
                        
                        for (int j = 0; j < maxStepsQ; ++j) {
                            const index_type indexQ = shellIndexQ + j;
                            
                            if ((shellIndexP != shellIndexQ) || (indexP >= indexQ)) {
                                const double coef = (indexP != indexQ) ? 2.0 : 1.0;
                                //const double coef = 2.0;
                                const double P_pq = coef * P.get(indexP, indexQ);
                                
                                for (int k = 0; k < maxStepsR; ++k) {
                                    const index_type indexR = shellIndexR + k;
                                    
                                    const double value = this->pEriEngines_[threadID].WORK[index];
                                    time_store.start();
                                    pRho->add(indexR, P_pq * value);
                                    time_store.stop();
                                    ++index;
                                }
                            } else {
                                index += maxStepsR;
                            }
                        }
                    }
                }
            }
        }

#ifdef _OPENMP        
        {
            const int numOfThreads = omp_get_num_threads();
            for (int thread = 0; thread < numOfThreads; ++thread) {
                if (thread == threadID) {
                    this->elapsetime_store_ += time_store.getElapseTime();
                }
#pragma omp barrier                
            }
        }
#else
        {
            this->elapsetime_store_ = time_store.getElapseTime();
        }
#endif // _OPENMP
    }
}


void DfEriX::getJ(const TlVector& rho, TlSymmetricMatrix* pJ)
{
    assert(pJ != NULL);
    this->elapsetime_store_ = 0.0;
    TlTime time_all;
    time_all.start();

    const TlOrbitalInfo orbitalInfo((*(this->pPdfParam_))["coordinates"],
                                    (*(this->pPdfParam_))["basis_sets"]);
    const TlOrbitalInfo_Density orbitalInfo_Density((*(this->pPdfParam_))["coordinates"],
                                                    (*(this->pPdfParam_))["basis_sets_j"]);
    const ShellArrayTable shellArrayTable = this->makeShellArrayTable(orbitalInfo);
    const ShellArrayTable shellArrayTable_Density = this->makeShellArrayTable(orbitalInfo_Density);

    //const TlSparseSymmetricMatrix schwarzTable = this->makeSchwarzTable(orbitalInfo);

    pJ->resize(this->m_nNumOfAOs);
    //pJ->zeroClear();
    
    this->createEngines();
    DfTaskCtrl* pDfTaskCtrl = this->getDfTaskCtrlObject();
    
    std::vector<DfTaskCtrl::Task2> taskList;
    bool hasTask = pDfTaskCtrl->getQueue2(orbitalInfo,
                                          true,
                                          this->grainSize_, &taskList, true);
    while (hasTask == true) {
        this->getJ_part(orbitalInfo,
                        orbitalInfo_Density,
                        shellArrayTable_Density,
                        taskList,
                        rho, pJ);

        hasTask = pDfTaskCtrl->getQueue2(orbitalInfo,
                                         true,
                                         this->grainSize_, &taskList);
    }

    this->finalize(pJ);

    delete pDfTaskCtrl;
    pDfTaskCtrl = NULL;
    const double time_ERI_all = this->destroyEngines();

    // statics report
    time_all.stop();
    {
        this->log_.info(TlUtils::format("all time: %16.1f sec.", time_all.getElapseTime()));
        this->log_.info(TlUtils::format(" store:    %16.1f sec.", this->elapsetime_store_));
        this->log_.info(TlUtils::format(" ERI time: %16.1f sec.", time_ERI_all));
    }
}


void DfEriX::getJ_part(const TlOrbitalInfo& orbitalInfo,
                       const TlOrbitalInfo_Density& orbitalInfo_Density,
                       const ShellArrayTable& shellArrayTable_Density,
                       const std::vector<DfTaskCtrl::Task2>& taskList,
                       const TlVector& rho, TlMatrixObject* pP)
{
    const int maxShellType = orbitalInfo.getMaxShellType();
    const int taskListSize = taskList.size();
    const double pairwisePGTO_cutoffThreshold = this->cutoffEpsilon3_;
    
#pragma omp parallel
    {
        TlTime time_store;
        int threadID = 0;
#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif // _OPENMP

        this->pEriEngines_[threadID].setPrimitiveLevelThreshold(this->cutoffEpsilon3_);

#pragma omp for schedule(runtime)
        for (int i = 0; i < taskListSize; ++i) {
            const index_type shellIndexP = taskList[i].shellIndex1;
            const index_type shellIndexQ = taskList[i].shellIndex2;
            const int shellTypeP = orbitalInfo.getShellType(shellIndexP);
            const int shellTypeQ = orbitalInfo.getShellType(shellIndexQ);
            
            const int maxStepsP = 2 * shellTypeP + 1;
            const int maxStepsQ = 2 * shellTypeQ + 1;
            const TlPosition posP = orbitalInfo.getPosition(shellIndexP);
            const TlPosition posQ = orbitalInfo.getPosition(shellIndexQ);
            const DfEriEngine::CGTO_Pair PQ = this->pEriEngines_[threadID].getCGTO_pair(orbitalInfo,
                                                                                        shellIndexP,
                                                                                        shellIndexQ,
                                                                                        pairwisePGTO_cutoffThreshold);
            const DfEriEngine::Query queryPQ(0, 0, shellTypeP, shellTypeQ);

            for (int shellTypeR = maxShellType -1; shellTypeR >= 0; --shellTypeR) {
                const int maxStepsR = 2 * shellTypeR + 1;
            
                const int shellTypeS = 0;
                const DfEriEngine::Query queryRS(0, 0, shellTypeR, shellTypeS);
                
                const std::size_t numOfShellArrayR = shellArrayTable_Density[shellTypeR].size();
                for (std::size_t indexR = 0; indexR < numOfShellArrayR; ++indexR) {
                    const index_type shellIndexR = shellArrayTable_Density[shellTypeR][indexR];
                    
                    const DfEriEngine::CGTO_Pair RS = 
                        this->pEriEngines_[threadID].getCGTO_pair(orbitalInfo_Density,
                                                                  shellIndexR,
                                                                  -1,
                                                                  pairwisePGTO_cutoffThreshold);
                    this->pEriEngines_[threadID].calc(queryPQ, queryRS, PQ, RS);
                
                    int index = 0;
                    for (int i = 0; i < maxStepsP; ++i) {
                        const index_type indexP = shellIndexP + i;
                        
                        for (int j = 0; j < maxStepsQ; ++j) {
                            const index_type indexQ = shellIndexQ + j;
                            
                            if ((shellIndexP != shellIndexQ) || (indexP >= indexQ)) {
                                double value = 0.0;
                                for (int k = 0; k < maxStepsR; ++k) {
                                    const index_type indexR = shellIndexR + k;
                                    
                                    value += rho.get(indexR) * this->pEriEngines_[threadID].WORK[index];
                                    ++index;
                                }

                                time_store.start();
                                pP->add(indexP, indexQ, value);
                                time_store.stop();
                            } else {
                                index += maxStepsR;
                            }
                        }
                    }
                }
            }
        }

#ifdef _OPENMP        
        {
            const int numOfThreads = omp_get_num_threads();
            for (int thread = 0; thread < numOfThreads; ++thread) {
                if (thread == threadID) {
                    this->elapsetime_store_ += time_store.getElapseTime();
                }
#pragma omp barrier
            }
        }
#else
        {
            this->elapsetime_store_ = time_store.getElapseTime();
        }
#endif // _OPENMP
    }
}


// TODO: OpenMP化
void DfEriX::getJpq(const TlSymmetricMatrix& P, TlSymmetricMatrix* pJ)
{
    assert(pJ != NULL);
    //this->clearCutoffStats();
    
    // カットオフ値の設定
    const double maxDeltaP = P.getMaxAbsoluteElement();
    if ((maxDeltaP > 0.0) && (maxDeltaP < 1.0)) {
        this->cutoffThreshold_ /= maxDeltaP;
        this->log_.info(TlUtils::format(" new cutoff threshold = % e\n",
                                        this->cutoffThreshold_));
    }

    // 本計算
    if (this->isDebugExactJ_ == true) {
        this->log_.info("calculate J using DEBUG engine.");
        this->getJpq_exact(P, pJ);
    } else {
        this->getJpq_integralDriven(P, pJ);
    }

    //this->cutoffReport();
}


void DfEriX::getJpq_exact(const TlSymmetricMatrix& P, TlSymmetricMatrix* pJ)
{
    assert(pJ != NULL);
    // const index_type numOfAOs = this->m_nNumOfAOs;

    assert(pJ != NULL);
    pJ->resize(this->m_nNumOfAOs);
    
    DfEriEngine engine;
    engine.setPrimitiveLevelThreshold(0.0);
    
    const TlOrbitalInfo orbitalInfo((*(this->pPdfParam_))["coordinates"],
                                    (*(this->pPdfParam_))["basis_sets"]);
    const ShellArrayTable shellArrayTable = this->makeShellArrayTable(orbitalInfo);
    const ShellPairArrayTable shellPairArrayTable = this->getShellPairArrayTable(shellArrayTable);
    
    for (int shellTypeP = DfEriX::MAX_SHELL_TYPE -1; shellTypeP >= 0; --shellTypeP) {
        const int maxStepsP = 2 * shellTypeP + 1;
        const ShellArray shellArrayP = shellArrayTable[shellTypeP];
        ShellArray::const_iterator pItEnd = shellArrayP.end();

        for (int shellTypeQ = DfEriX::MAX_SHELL_TYPE -1; shellTypeQ >= 0; --shellTypeQ) {
            const int maxStepsQ = 2 * shellTypeQ + 1;
            const ShellArray shellArrayQ = shellArrayTable[shellTypeQ];
            ShellArray::const_iterator qItEnd = shellArrayQ.end();

            const DfEriEngine::Query queryPQ(0, 0, shellTypeP, shellTypeQ);

            for (int shellTypeR = DfEriX::MAX_SHELL_TYPE -1; shellTypeR >= 0; --shellTypeR) {
                const int maxStepsR = 2 * shellTypeR + 1;
                const ShellArray shellArrayR = shellArrayTable[shellTypeR];
                ShellArray::const_iterator rItEnd = shellArrayR.end();

                for (int shellTypeS = DfEriX::MAX_SHELL_TYPE -1; shellTypeS >= 0; --shellTypeS) {
                    const int maxStepsS = 2 * shellTypeS + 1;
                    const ShellArray shellArrayS = shellArrayTable[shellTypeS];
                    ShellArray::const_iterator sItEnd = shellArrayS.end();

                    const DfEriEngine::Query queryRS(0, 0, shellTypeR, shellTypeS);
        
                    for (ShellArray::const_iterator pIt = shellArrayP.begin(); pIt != pItEnd; ++pIt) {
                        const index_type shellIndexP = *pIt;
                        for (ShellArray::const_iterator qIt = shellArrayQ.begin(); qIt != qItEnd; ++qIt) {
                            const index_type shellIndexQ = *qIt;

                            const DfEriEngine::CGTO_Pair PQ = engine.getCGTO_pair(orbitalInfo,
                                                                                  shellIndexP, shellIndexQ,
                                                                                  0.0);
                            
                            for (ShellArray::const_iterator rIt = shellArrayR.begin(); rIt != rItEnd; ++rIt) {
                                const index_type shellIndexR = *rIt;
                                for (ShellArray::const_iterator sIt = shellArrayS.begin(); sIt != sItEnd; ++sIt) {
                                    const index_type shellIndexS = *sIt;

                                    const DfEriEngine::CGTO_Pair RS = engine.getCGTO_pair(orbitalInfo,
                                                                                          shellIndexR, shellIndexS,
                                                                                          0.0);
                                    engine.calc(queryPQ, queryRS, PQ, RS);

                                    int index = 0;
                                    for (int i = 0; i < maxStepsP; ++i) {
                                        const int indexP = shellIndexP + i;
                                        for (int j = 0; j < maxStepsQ; ++j) {
                                            const int indexQ = shellIndexQ + j;
                                            
                                            for (int k = 0; k < maxStepsR; ++k) {
                                                const int indexR = shellIndexR + k;
                                                for (int l = 0; l < maxStepsS; ++l) {
                                                    const int indexS = shellIndexS + l;
                                                    
                                                    if (indexP >= indexQ) {
                                                        const double P_rs = P.get(indexR, indexS);
                                                        const double value = engine.WORK[index];
                                                        pJ->add(indexP, indexQ, P_rs * value);
                                                    }
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
}


void DfEriX::getJpq_integralDriven(const TlSymmetricMatrix& P, TlSymmetricMatrix* pJ)
{
    assert(pJ != NULL);
    this->elapsetime_store_ = 0.0;
    TlTime time_all;
    time_all.start();

    pJ->resize(this->m_nNumOfAOs);

    const TlOrbitalInfo orbitalInfo((*(this->pPdfParam_))["coordinates"],
                                    (*(this->pPdfParam_))["basis_sets"]);

    const TlSparseSymmetricMatrix schwarzTable = this->makeSchwarzTable(orbitalInfo);

#ifdef DEBUG_J
    const index_type numOfAOs = this->m_nNumOfAOs;
    this->IA_J_ID1_.resize(numOfAOs);
    this->IA_J_ID2_.resize(numOfAOs);
#endif // DEBUG_J
    
    this->createEngines();
    DfTaskCtrl* pDfTaskCtrl = this->getDfTaskCtrlObject();
    pDfTaskCtrl->setCutoffThreshold(this->cutoffThreshold_);
    pDfTaskCtrl->setCutoffEpsilon_density(0.0);  // cannot use this cutoff
    pDfTaskCtrl->setCutoffEpsilon_distribution(this->cutoffEpsilon_distribution_);

    std::vector<DfTaskCtrl::Task4> taskList;
    bool hasTask = pDfTaskCtrl->getQueue4(orbitalInfo,
                                          schwarzTable,
                                          this->grainSize_,
                                          &taskList,
                                          true);
    while (hasTask == true) {
        this->getJ_integralDriven_part(orbitalInfo,
                                       taskList,
                                       P, pJ);
        hasTask = pDfTaskCtrl->getQueue4(orbitalInfo,
                                         schwarzTable,
                                         this->grainSize_,
                                         &taskList);
    }
                    
    this->finalize(pJ);
    //pDfTaskCtrl->cutoffReport();

    delete pDfTaskCtrl;
    pDfTaskCtrl = NULL;
    const double time_ERI_all = this->destroyEngines();
    
    // statics report
    time_all.stop();
    {
        this->log_.info(TlUtils::format("all time:  %16.1f sec.", time_all.getElapseTime()));
        this->log_.info(TlUtils::format(" store:    %16.1f sec.", this->elapsetime_store_));
        this->log_.info(TlUtils::format(" ERI time: %16.1f sec.", time_ERI_all));
    }

    // debug
#ifdef DEBUG_J
    if (this->isDebugOutJ_ == true) {
        for (index_type i = 0; i < numOfAOs; ++i) {
            for (index_type j = 0; j <= i; ++j) {
                
                std::cerr << TlUtils::format(">>>>J(%2d,%2d)", i, j)
                          << std::endl;
                for (index_type k = 0; k < numOfAOs; ++k) {
                    for (index_type l = 0; l <= k; ++l) {
                        const int counter1 = this->IA_J_ID1_.getCount(i, j, k, l);
                        const int counter2 = this->IA_J_ID2_.getCount(i, j, k, l);
                        const int counter = counter1 + counter2;
                        std::string YN = "  ";
                        if (counter != ((k == l) ? 1 : 2)) {
                            YN = "NG";
                        }
                        std::cerr << TlUtils::format("J(%2d,%2d) <= (%2d,%2d) %2d(%2d,%2d) %s",
                                                     i, j, k, l, counter, counter1, counter2, YN.c_str())
                                  << std::endl;
                    }
                }
                std::cerr << std::endl;
            }
        }
    }
#endif // DEBUG_J
}


void DfEriX::getJ_integralDriven_part(const TlOrbitalInfoObject& orbitalInfo,
                                      const std::vector<DfTaskCtrl::Task4>& taskList,
                                      const TlMatrixObject& P, TlMatrixObject* pJ)
{
    const int taskListSize = taskList.size();
    const double pairwisePGTO_cutoffThreshold = this->cutoffEpsilon3_;

#pragma omp parallel
    {
        TlTime time_store;
        int threadID = 0;

#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif // _OPENMP
        
        this->pEriEngines_[threadID].setPrimitiveLevelThreshold(this->cutoffEpsilon3_);

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
                        
            const DfEriEngine::CGTO_Pair PQ = this->pEriEngines_[threadID].getCGTO_pair(orbitalInfo,
                                                                                        shellIndexP,
                                                                                        shellIndexQ,
                                                                                        pairwisePGTO_cutoffThreshold);
            const DfEriEngine::CGTO_Pair RS = this->pEriEngines_[threadID].getCGTO_pair(orbitalInfo,
                                                                                        shellIndexR,
                                                                                        shellIndexS,
                                                                                        pairwisePGTO_cutoffThreshold);
            const DfEriEngine::Query queryPQ(0, 0, shellTypeP, shellTypeQ);
            const DfEriEngine::Query queryRS(0, 0, shellTypeR, shellTypeS);

            this->pEriEngines_[threadID].calc(queryPQ, queryRS, PQ, RS);
                        
            this->storeJ_integralDriven(shellIndexP, maxStepsP,
                                        shellIndexQ, maxStepsQ,
                                        shellIndexR, maxStepsR,
                                        shellIndexS, maxStepsS,
                                        this->pEriEngines_[threadID], P, pJ);
        }

#ifdef _OPENMP        
        {
            const int numOfThreads = omp_get_num_threads();
            for (int thread = 0; thread < numOfThreads; ++thread) {
                if (thread == threadID) {
                    this->elapsetime_store_ += time_store.getElapseTime();
                }
#pragma omp barrier
            }
        }
#else
        {
            this->elapsetime_store_ = time_store.getElapseTime();
        }
#endif // _OPENMP
    }
        
}


void DfEriX::storeJ_integralDriven(const index_type shellIndexP, const int maxStepsP,
                                   const index_type shellIndexQ, const int maxStepsQ,
                                   const index_type shellIndexR, const int maxStepsR,
                                   const index_type shellIndexS, const int maxStepsS,
                                   const DfEriEngine& engine,
                                   const TlMatrixObject& P,
                                   TlMatrixObject* pJ)
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
                        pJ->add(indexP, indexQ, coefEq1 * P_rs * value);
#ifdef DEBUG_J
                        this->IA_J_ID1_.countUp(indexP, indexQ, indexR, indexS, coefEq1);
#endif // DEBUG_J
                        
                        // Eq.2 : (indexR, indexS) <= (indexP, indexQ)
                        if ((shellIndexP != shellIndexR) || (shellIndexQ != shellIndexS) || (indexP == indexR)) {
                            if (((indexP + indexQ) != (indexR + indexS)) ||
                                ((indexP * indexQ) != (indexR * indexS))) {
                                // Eq.1の条件と重複しないようにするための措置
                                            
                                const double coefEq2 = (indexP != indexQ) ? 2.0 : 1.0;
                                pJ->add(indexR, indexS, coefEq2 * P_pq * value);
#ifdef DEBUG_J
                                this->IA_J_ID2_.countUp(indexR, indexS, indexP, indexQ, coefEq2);
#endif // DEBUG_J
                            }
                        }
                    }
                    ++index;
                }
            }
        }
    }
}


void DfEriX::getJab(TlSymmetricMatrix* pJab)
{
    assert(pJab != NULL);
    const index_type numOfAuxDens = this->m_nNumOfAux;
    pJab->resize(numOfAuxDens);

    const TlOrbitalInfo_Density orbitalInfo_Density((*(this->pPdfParam_))["coordinates"],
                                                    (*(this->pPdfParam_))["basis_sets_j"]);

    const ShellArrayTable shellArrayTable = this->makeShellArrayTable(orbitalInfo_Density);

    this->createEngines();
    DfTaskCtrl* pDfTaskCtrl = this->getDfTaskCtrlObject();

    std::vector<DfTaskCtrl::Task2> taskList;
    pDfTaskCtrl->setCutoffThreshold(this->cutoffThreshold_);
    pDfTaskCtrl->setCutoffEpsilon_density(0.0);  // cannot use this cutoff
    pDfTaskCtrl->setCutoffEpsilon_distribution(this->cutoffEpsilon_distribution_);

    bool hasTask = pDfTaskCtrl->getQueue2(orbitalInfo_Density,
                                          false,
                                          this->grainSize_, &taskList, true);
    while (hasTask == true) {
        this->getJab_part(orbitalInfo_Density,
                          taskList,
                          pJab);

        hasTask = pDfTaskCtrl->getQueue2(orbitalInfo_Density,
                                         false,
                                         this->grainSize_, &taskList);
    }

    this->finalize(pJab);

    delete pDfTaskCtrl;
    pDfTaskCtrl = NULL;
    this->destroyEngines();
}


void DfEriX::getJab_part(const TlOrbitalInfoObject& orbitalInfo,
                         const std::vector<DfTaskCtrl::Task2>& taskList,
                         TlMatrixObject* pJab)
{
    //const int maxShellType = orbitalInfo.getMaxShellType();
    const int taskListSize = taskList.size();
    const double pairwisePGTO_cutoffThreshold = this->cutoffEpsilon3_;

#pragma omp parallel
    {
        int threadID = 0;
#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif // _OPENMP

        this->pEriEngines_[threadID].setPrimitiveLevelThreshold(this->cutoffEpsilon3_);
        
#pragma omp for schedule(runtime)
        for (int i = 0; i < taskListSize; ++i) {
            const index_type shellIndexP = taskList[i].shellIndex1;
            const index_type shellIndexR = taskList[i].shellIndex2;
            const int shellTypeP = orbitalInfo.getShellType(shellIndexP);
            const int shellTypeR = orbitalInfo.getShellType(shellIndexR);

            // schwarz cutoff
            // const int shellQuartetType =
            //     ((shellTypeP * maxShellType + shellTypeQ) * maxShellType + shellTypeP) * maxShellType + shellTypeQ;
            // const bool isAlive = this->isAliveBySchwarzCutoff(shellIndexP, shellIndexQ,
            //                                                   shellIndexP, shellIndexQ,
            //                                                   shellQuartetType,
            //                                                   schwarzTable,
            //                                                   this->cutoffThreshold_);
            // if (isAlive != true) {
            //     continue;
            // }

            const int maxStepsP = 2 * shellTypeP + 1;
            const int maxStepsR = 2 * shellTypeR + 1;
            //const TlPosition posP = orbitalInfo.getPosition(shellIndexP);
            //const TlPosition posR = orbitalInfo.getPosition(shellIndexR);
            const DfEriEngine::Query queryPQ(0, 0, shellTypeP, 0);
            const DfEriEngine::Query queryRS(0, 0, shellTypeR, 0);
            const DfEriEngine::CGTO_Pair PQ = this->pEriEngines_[threadID].getCGTO_pair(orbitalInfo,
                                                                                        shellIndexP,
                                                                                        -1,
                                                                                        pairwisePGTO_cutoffThreshold);
            const DfEriEngine::CGTO_Pair RS = this->pEriEngines_[threadID].getCGTO_pair(orbitalInfo,
                                                                                        shellIndexR,
                                                                                        -1,
                                                                                        pairwisePGTO_cutoffThreshold);

            this->pEriEngines_[threadID].calc(queryPQ, queryRS, PQ, RS);
                
            int index = 0;
            for (int p = 0; p < maxStepsP; ++p) {
                const index_type indexP = shellIndexP + p;
                
                for (int r = 0; r < maxStepsR; ++r) {
                    const index_type indexR = shellIndexR + r;

                    if ((shellIndexP != shellIndexR) || (indexP >= indexR)) {
                        const double value = this->pEriEngines_[threadID].WORK[index];
                        pJab->add(indexP, indexR, value);
                    }
                    ++index;
                }
            }
        }
    }
}


void DfEriX::getForceJ(const TlSymmetricMatrix& P, TlMatrix* pForce)
{
    assert(pForce != NULL);
    pForce->resize(this->m_nNumOfAtoms, 3);
    
    const TlOrbitalInfo orbitalInfo((*(this->pPdfParam_))["coordinates"],
                                    (*(this->pPdfParam_))["basis_sets"]);
    const TlSparseSymmetricMatrix schwarzTable = this->makeSchwarzTable(orbitalInfo);

    this->createEngines();
    DfTaskCtrl* pDfTaskCtrl = this->getDfTaskCtrlObject();

    std::vector<DfTaskCtrl::Task4> taskList;
    pDfTaskCtrl->setCutoffThreshold(this->cutoffThreshold_);
    pDfTaskCtrl->setCutoffEpsilon_density(0.0);  // cannot use this cutoff
    pDfTaskCtrl->setCutoffEpsilon_distribution(this->cutoffEpsilon_distribution_);

    bool hasTask = pDfTaskCtrl->getQueue_Force4(orbitalInfo,
                                                schwarzTable,
                                                this->grainSize_, &taskList, true);
    while (hasTask == true) {
        this->getForceJ_part(orbitalInfo, taskList,
                             P, pForce);
        hasTask = pDfTaskCtrl->getQueue_Force4(orbitalInfo,
                                               schwarzTable,
                                               this->grainSize_, &taskList);
    }

    this->finalize(pForce);

    delete pDfTaskCtrl;
    pDfTaskCtrl = NULL;
    this->destroyEngines();
}


void DfEriX::getForceJ_part(const TlOrbitalInfoObject& orbitalInfo,
                            const std::vector<DfTaskCtrl::Task4>& taskList,
                            const TlMatrixObject& P, TlMatrixObject* pForce)
{
    const int taskListSize = taskList.size();
    const double pairwisePGTO_cutoffThreshold = this->cutoffEpsilon3_;
    
#pragma omp parallel
    {
        int threadID = 0;
#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif // _OPENMP
        
        this->pEriEngines_[threadID].setPrimitiveLevelThreshold(this->cutoffEpsilon3_);

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
            const index_type atomIndexA = orbitalInfo.getAtomIndex(shellIndexP);
            const index_type atomIndexB = orbitalInfo.getAtomIndex(shellIndexQ);
            const index_type atomIndexC = orbitalInfo.getAtomIndex(shellIndexR);
            const index_type atomIndexD = orbitalInfo.getAtomIndex(shellIndexS);

            // std::cerr << TlUtils::format("[%d %d|%d %d]",
            //                              shellIndexP, shellIndexQ, shellIndexR, shellIndexS)
            //           << std::endl;
            
            if ((atomIndexA == atomIndexB) && (atomIndexB == atomIndexC) &&
                (atomIndexC == atomIndexD) && (atomIndexD == atomIndexA)) {
                continue;
            }
            
            const DfEriEngine::Query queryPQ00(0, 0, shellTypeP, shellTypeQ);
            const DfEriEngine::Query queryRS00(0, 0, shellTypeR, shellTypeS);
            const DfEriEngine::CGTO_Pair PQ =
                this->pEriEngines_[threadID].getCGTO_pair(orbitalInfo,
                                                          shellIndexP,
                                                          shellIndexQ,
                                                          pairwisePGTO_cutoffThreshold);
            const DfEriEngine::CGTO_Pair RS =
                this->pEriEngines_[threadID].getCGTO_pair(orbitalInfo,
                                                          shellIndexR,
                                                          shellIndexS,
                                                          pairwisePGTO_cutoffThreshold);
            
            this->pEriEngines_[threadID].calcGrad(queryPQ00, queryRS00, PQ, RS);
                        
            this->storeForceJ_integralDriven(atomIndexA, atomIndexB,
                                             atomIndexC, atomIndexD,
                                             shellIndexP, maxStepsP,
                                             shellIndexQ, maxStepsQ,
                                             shellIndexR, maxStepsR,
                                             shellIndexS, maxStepsS,
                                             this->pEriEngines_[threadID], P, pForce);
        }
    }
}


void DfEriX::storeForceJ_integralDriven(const int atomIndexA, const int atomIndexB,
                                        const int atomIndexC, const int atomIndexD,
                                        const index_type shellIndexP, const int maxStepsP,
                                        const index_type shellIndexQ, const int maxStepsQ,
                                        const index_type shellIndexR, const int maxStepsR,
                                        const index_type shellIndexS, const int maxStepsS,
                                        const DfEriEngine& engine,
                                        const TlMatrixObject& P,
                                        TlMatrixObject* pForce)
{
    int index = 0;
    this->storeForceJ_integralDriven(atomIndexA, atomIndexB,
                                     atomIndexC, atomIndexD,
                                     shellIndexP, maxStepsP,
                                     shellIndexQ, maxStepsQ,
                                     shellIndexR, maxStepsR,
                                     shellIndexS, maxStepsS,
                                     engine, P, pForce,
                                     X, &index);
    this->storeForceJ_integralDriven(atomIndexA, atomIndexB,
                                     atomIndexC, atomIndexD,
                                     shellIndexP, maxStepsP,
                                     shellIndexQ, maxStepsQ,
                                     shellIndexR, maxStepsR,
                                     shellIndexS, maxStepsS,
                                     engine, P, pForce,
                                     Y, &index);
    this->storeForceJ_integralDriven(atomIndexA, atomIndexB,
                                     atomIndexC, atomIndexD,
                                     shellIndexP, maxStepsP,
                                     shellIndexQ, maxStepsQ,
                                     shellIndexR, maxStepsR,
                                     shellIndexS, maxStepsS,
                                     engine, P, pForce,
                                     Z, &index);
}


void DfEriX::storeForceJ_integralDriven(const int atomIndexA, const int atomIndexB,
                                        const int atomIndexC, const int atomIndexD,
                                        const index_type shellIndexP, const int maxStepsP,
                                        const index_type shellIndexQ, const int maxStepsQ,
                                        const index_type shellIndexR, const int maxStepsR,
                                        const index_type shellIndexS, const int maxStepsS,
                                        const DfEriEngine& engine,
                                        const TlMatrixObject& P,
                                        TlMatrixObject* pForce,
                                        const int target, int* pIndex)
{
    for (int stepP = 0; stepP < maxStepsP; ++stepP) {
        const index_type indexP = shellIndexP + stepP;
        const index_type iw = indexP * (indexP -1) / 2;

        for (int stepQ = 0; stepQ < maxStepsQ; ++stepQ) {
            const index_type indexQ = shellIndexQ + stepQ;

            if (indexQ <= indexP) {
                const double cij = (indexP != indexQ) ? 2.0 : 1.0;
                const index_type ij = iw + indexQ;
                const double dcij = P.get(indexP, indexQ);
                
                for (int stepR = 0; stepR < maxStepsR; ++stepR) {
                    const index_type indexR = shellIndexR + stepR;
                    const index_type kw = indexR * (indexR -1) / 2;
                    
                    for (int stepS = 0; stepS < maxStepsS; ++stepS) {
                        const index_type indexS = shellIndexS + stepS;
                        const index_type kl = kw + indexS;
                        
                        if ((indexS <= indexR) && (ij >= kl)) {
                            double cijkl = cij * ((indexR != indexS) ? 2.0 : 1.0);
                            cijkl *= (ij != kl) ? 2.0 : 1.0;

                            const double coef = cijkl * dcij * P.get(indexR, indexS);
                            const double gradIntA = engine.WORK_A[*pIndex];
                            const double gradIntB = engine.WORK_B[*pIndex];
                            const double gradIntC = engine.WORK_C[*pIndex];
                            const double gradIntD = - (gradIntA + gradIntB + gradIntC);
                            
                            pForce->add(atomIndexA, target, coef * gradIntA);
                            pForce->add(atomIndexB, target, coef * gradIntB);
                            pForce->add(atomIndexC, target, coef * gradIntC);
                            pForce->add(atomIndexD, target, coef * gradIntD);
                        }
                        
                        ++(*pIndex);
                    }
                }
            } else {
                *pIndex += (maxStepsR * maxStepsS);
            }
        }
    }
}


void DfEriX::getForceJ(const TlSymmetricMatrix& P, const TlVector& rho,
                       TlMatrix* pForce)
{
    assert(pForce != NULL);
    //this->clearCutoffStats();

    // カットオフ値の設定
    const double maxDeltaP = P.getMaxAbsoluteElement();
    if (maxDeltaP < 1.0) {
        this->cutoffThreshold_ /= std::fabs(maxDeltaP);
        this->log_.info(TlUtils::format(" new cutoff threshold = % e",
                                        this->cutoffThreshold_));
    }
    
    pForce->resize(this->m_nNumOfAtoms, 3);
    
    const TlOrbitalInfo orbitalInfo((*(this->pPdfParam_))["coordinates"],
                                    (*(this->pPdfParam_))["basis_sets"]);
    const TlOrbitalInfo_Density orbitalInfo_Density((*(this->pPdfParam_))["coordinates"],
                                                    (*(this->pPdfParam_))["basis_sets_j"]);
    const ShellArrayTable shellArrayTable = this->makeShellArrayTable(orbitalInfo);
    const ShellArrayTable shellArrayTable_Density = this->makeShellArrayTable(orbitalInfo_Density);

    this->createEngines();
    DfTaskCtrl* pDfTaskCtrl = this->getDfTaskCtrlObject();
    pDfTaskCtrl->setCutoffThreshold(this->cutoffThreshold_);
    pDfTaskCtrl->setCutoffEpsilon_density(0.0);  // cannot use this cutoff
    pDfTaskCtrl->setCutoffEpsilon_distribution(this->cutoffEpsilon_distribution_);

    std::vector<DfTaskCtrl::Task2> taskList;
    bool hasTask = pDfTaskCtrl->getQueue2(orbitalInfo,
                                          false,
                                          this->grainSize_, &taskList, true);
    while (hasTask == true) {
        this->getForceJ_part(orbitalInfo,
                             orbitalInfo_Density,
                             shellArrayTable_Density,
                             taskList,
                             P, rho, pForce);

        hasTask = pDfTaskCtrl->getQueue2(orbitalInfo,
                                         false,
                                         this->grainSize_, &taskList);
    }

    this->finalize(pForce);

    delete pDfTaskCtrl;
    pDfTaskCtrl = NULL;
    this->destroyEngines();
}


void DfEriX::getForceJ_part(const TlOrbitalInfoObject& orbitalInfo,
                            const TlOrbitalInfoObject& orbitalInfo_Density,
                            const ShellArrayTable& shellArrayTable_Density,
                            std::vector<DfTaskCtrl::Task2>& taskList,
                            const TlSymmetricMatrix& P, const TlVector& rho,
                            TlMatrix* pForce)
{
    static const int BUFFER_SIZE = 3 * 5 * 5 * 5; // (xyz) * 5d * 5d * 5d
    const int maxShellType = orbitalInfo.getMaxShellType();
    const int taskListSize = taskList.size();
    const double pairwisePGTO_cutoffThreshold = this->cutoffEpsilon3_;

#pragma omp parallel
    {
        int threadID = 0;
#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif // _OPENMP

        this->pEriEngines_[threadID].setPrimitiveLevelThreshold(this->cutoffEpsilon3_);

        double* p_dJdA = new double[BUFFER_SIZE];
        double* p_dJdB = new double[BUFFER_SIZE];

#pragma omp for schedule(runtime)
        for (int i = 0; i < taskListSize; ++i) {
            const index_type shellIndexP = taskList[i].shellIndex1;
            const index_type shellIndexQ = taskList[i].shellIndex2;
            const int shellTypeP = orbitalInfo.getShellType(shellIndexP);
            const int shellTypeQ = orbitalInfo.getShellType(shellIndexQ);
            const int maxStepsP = 2 * shellTypeP + 1;
            const int maxStepsQ = 2 * shellTypeQ + 1;
            const TlPosition posP = orbitalInfo.getPosition(shellIndexP);
            const TlPosition posQ = orbitalInfo.getPosition(shellIndexQ);
            const index_type atomIndexA = orbitalInfo.getAtomIndex(shellIndexP);
            const index_type atomIndexB = orbitalInfo.getAtomIndex(shellIndexQ);
            const DfEriEngine::CGTO_Pair PQ = this->pEriEngines_[threadID].getCGTO_pair(orbitalInfo,
                                                                                        shellIndexP,
                                                                                        shellIndexQ,
                                                                                        pairwisePGTO_cutoffThreshold);
            const DfEriEngine::Query queryPQ(1, 0, shellTypeP, shellTypeQ);
            const DfEriEngine::Query queryQP(0, 1, shellTypeP, shellTypeQ);

            for (int shellTypeR = maxShellType -1; shellTypeR >= 0; --shellTypeR) {
                const int maxStepsR = 2 * shellTypeR + 1;
                const ShellArray& shellArrayR = shellArrayTable_Density[shellTypeR];
                const DfEriEngine::Query queryRS(0, 0, shellTypeR, 0);

                const std::size_t shellArraySizeR = shellArrayR.size();
                for (std::size_t r = 0; r < shellArraySizeR; ++r) {
                    const index_type shellIndexR = shellArrayR[r];
                    const index_type atomIndexC = orbitalInfo_Density.getAtomIndex(shellIndexR);

                    if ((atomIndexA == atomIndexB) && (atomIndexB == atomIndexC) &&
                        (atomIndexC == atomIndexA)) {
                        continue;
                    }
                    
                    const DfEriEngine::CGTO_Pair RS = this->pEriEngines_[threadID].getCGTO_pair(orbitalInfo_Density,
                                                                                                shellIndexR,
                                                                                                -1,
                                                                                                pairwisePGTO_cutoffThreshold);
                    this->pEriEngines_[threadID].calc(queryPQ, queryRS, PQ, RS);
                    std::copy(this->pEriEngines_[threadID].WORK,
                              this->pEriEngines_[threadID].WORK + BUFFER_SIZE, p_dJdA);
                    this->pEriEngines_[threadID].calc(queryQP, queryRS, PQ, RS);
                    std::copy(this->pEriEngines_[threadID].WORK,
                              this->pEriEngines_[threadID].WORK + BUFFER_SIZE, p_dJdB);
                    
                    int index = 0;
                    this->storeForceJ(atomIndexA, atomIndexB, atomIndexC,
                                      shellIndexP, maxStepsP,
                                      shellIndexQ, maxStepsQ,
                                      shellIndexR, maxStepsR,
                                      p_dJdA, p_dJdB, P, rho,
                                      pForce, X, &index);
                    this->storeForceJ(atomIndexA, atomIndexB, atomIndexC,
                                      shellIndexP, maxStepsP,
                                      shellIndexQ, maxStepsQ,
                                      shellIndexR, maxStepsR,
                                      p_dJdA, p_dJdB, P, rho,
                                      pForce, Y, &index);
                    this->storeForceJ(atomIndexA, atomIndexB, atomIndexC,
                                      shellIndexP, maxStepsP,
                                      shellIndexQ, maxStepsQ,
                                      shellIndexR, maxStepsR,
                                      p_dJdA, p_dJdB, P, rho,
                                      pForce, Z, &index);

                }
            }
        }

        delete[] p_dJdA;
        p_dJdA = NULL;
        delete[] p_dJdB;
        p_dJdB = NULL;
    }
}


void DfEriX::storeForceJ(const index_type atomIndexA,
                         const index_type atomIndexB,
                         const index_type atomIndexC,
                         const index_type shellIndexP, const int maxStepsP,
                         const index_type shellIndexQ, const int maxStepsQ,
                         const index_type shellIndexR, const int maxStepsR,
                         const double* p_dJdA, const double* p_dJdB,
                         const TlMatrixObject& P,
                         const TlVectorObject& rho,
                         TlMatrixObject* pForce,
                         const int target, int* pIndex)
{
    for (int stepP = 0; stepP < maxStepsP; ++stepP) {
        const index_type indexP = shellIndexP + stepP;

        for (int stepQ = 0; stepQ < maxStepsQ; ++stepQ) {
            const index_type indexQ = shellIndexQ + stepQ;

            if ((shellIndexP != shellIndexQ) || (indexP >= indexQ)) {
                double coef_PQ = P.get(indexP, indexQ);
                coef_PQ *= (indexP != indexQ) ? 2.0 : 1.0;

                for (int stepR = 0; stepR < maxStepsR; ++stepR) {
                    const double gradIntA = p_dJdA[*pIndex];
                    const double gradIntB = p_dJdB[*pIndex];
                    const double gradIntC = - (gradIntA + gradIntB);
                    const double coef_PQR = coef_PQ * rho[shellIndexR + stepR];

                    pForce->add(atomIndexA, target, coef_PQR * gradIntA);
                    pForce->add(atomIndexB, target, coef_PQR * gradIntB);
                    pForce->add(atomIndexC, target, coef_PQR * gradIntC);
                    ++(*pIndex);
                }
            } else {
                *pIndex += maxStepsR;
            }
        }
    }
}


void DfEriX::getForceJ(const TlVector& rho, TlMatrix* pForce)
{
    assert(pForce != NULL);
    pForce->resize(this->m_nNumOfAtoms, 3);
    
    const TlOrbitalInfo_Density orbitalInfo_Density((*(this->pPdfParam_))["coordinates"],
                                                    (*(this->pPdfParam_))["basis_sets_j"]);
    const ShellArrayTable shellArray_Density = this->makeShellArrayTable(orbitalInfo_Density);

    this->createEngines();
    DfTaskCtrl* pDfTaskCtrl = this->getDfTaskCtrlObject();
    pDfTaskCtrl->setCutoffThreshold(this->cutoffThreshold_);
    pDfTaskCtrl->setCutoffEpsilon_density(0.0);  // cannot use this cutoff
    pDfTaskCtrl->setCutoffEpsilon_distribution(this->cutoffEpsilon_distribution_);

    std::vector<DfTaskCtrl::Task2> taskList;
    bool hasTask = pDfTaskCtrl->getQueue2(orbitalInfo_Density,
                                          false,
                                          this->grainSize_, &taskList, true);
    while (hasTask == true) {
        this->getForceJ_part(orbitalInfo_Density,
                             taskList,
                             rho, pForce);

        hasTask = pDfTaskCtrl->getQueue2(orbitalInfo_Density,
                                         false,
                                         this->grainSize_, &taskList);
    }

    this->finalize(pForce);

    delete pDfTaskCtrl;
    pDfTaskCtrl = NULL;
    this->destroyEngines();
}


void DfEriX::getForceJ_part(const TlOrbitalInfoObject& orbitalInfo_Density,
                            std::vector<DfTaskCtrl::Task2>& taskList,
                            const TlVector& rho,
                            TlMatrix* pForce)
{
    //static const int BUFFER_SIZE = 3 * 5 * 5 * 5; // (xyz) * 5d * 5d * 5d
    //const int maxShellType = orbitalInfo.getMaxShellType();
    const int taskListSize = taskList.size();
    const double pairwisePGTO_cutoffThreshold = this->cutoffEpsilon3_;

#pragma omp parallel
    {
        int threadID = 0;
#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif // _OPENMP

        this->pEriEngines_[threadID].setPrimitiveLevelThreshold(this->cutoffEpsilon3_);

#pragma omp for schedule(runtime)
        for (int i = 0; i < taskListSize; ++i) {
            const index_type shellIndexP = taskList[i].shellIndex1;
            const index_type shellIndexR = taskList[i].shellIndex2;
            const index_type atomIndexA = orbitalInfo_Density.getAtomIndex(shellIndexP);
            const index_type atomIndexC = orbitalInfo_Density.getAtomIndex(shellIndexR);
            if (atomIndexA == atomIndexC) {
                continue;
            }

            const int shellTypeP = orbitalInfo_Density.getShellType(shellIndexP);
            const int shellTypeR = orbitalInfo_Density.getShellType(shellIndexR);
            const int maxStepsP = 2 * shellTypeP + 1;
            const int maxStepsR = 2 * shellTypeR + 1;
            const DfEriEngine::Query queryPQ(1, 0, shellTypeP, 0);
            const DfEriEngine::Query queryRS(0, 0, shellTypeR, 0);
            const DfEriEngine::CGTO_Pair PQ = this->pEriEngines_[threadID].getCGTO_pair(orbitalInfo_Density,
                                                                                        shellIndexP,
                                                                                        -1,
                                                                                        pairwisePGTO_cutoffThreshold);
            const DfEriEngine::CGTO_Pair RS = this->pEriEngines_[threadID].getCGTO_pair(orbitalInfo_Density,
                                                                                        shellIndexR,
                                                                                        -1,
                                                                                        pairwisePGTO_cutoffThreshold);

            this->pEriEngines_[threadID].calc(queryPQ, queryRS, PQ, RS);

            int index = 0;
            this->storeForceJ(atomIndexA, atomIndexC,
                              shellIndexP, maxStepsP,
                              shellIndexR, maxStepsR,
                              this->pEriEngines_[threadID],
                              rho, pForce,
                              X, &index);
            this->storeForceJ(atomIndexA, atomIndexC,
                              shellIndexP, maxStepsP,
                              shellIndexR, maxStepsR,
                              this->pEriEngines_[threadID],
                              rho, pForce,
                              Y, &index);
            this->storeForceJ(atomIndexA, atomIndexC,
                              shellIndexP, maxStepsP,
                              shellIndexR, maxStepsR,
                              this->pEriEngines_[threadID],
                              rho, pForce,
                              Z, &index);
        }
    }
}


void DfEriX::storeForceJ(const index_type atomIndexA,
                         const index_type atomIndexC,
                         const index_type shellIndexP, const int maxStepsP,
                         const index_type shellIndexR, const int maxStepsR,
                         const DfEriEngine& engine,
                         const TlVectorObject& rho,
                         TlMatrixObject* pForce,
                         const int target, int* pIndex)
{
    for (int stepP = 0; stepP < maxStepsP; ++stepP) {
        const index_type indexP = shellIndexP + stepP;
        const double coef_P = rho[indexP];
        
        for (int stepR = 0; stepR < maxStepsR; ++stepR) {
            const index_type indexR = shellIndexR + stepR;

            if ((shellIndexP != shellIndexR) || (indexP >= indexR)) {
                double coef = coef_P * rho[indexR];
                coef *= (indexP != indexR) ? 2.0 : 1.0;
                const double gradA = engine.WORK[*pIndex];
                const double gradC = - gradA;
                
                pForce->add(atomIndexA, target, coef * gradA);
                pForce->add(atomIndexC, target, coef * gradC);
            }
            ++(*pIndex);
        }
    }
}


void DfEriX::getK(const TlSymmetricMatrix& P, TlSymmetricMatrix* pK)
{
    //this->clearCutoffStats();

    // カットオフ値の設定
    const double maxDeltaP = P.getMaxAbsoluteElement();
    if ((maxDeltaP > 0.0) && (maxDeltaP < 1.0)) {
        this->cutoffThreshold_ /= maxDeltaP;
        this->log_.info(TlUtils::format(" new cutoff threshold = % e",
                                        this->cutoffThreshold_));
    }

    // 本計算
    if (this->isDebugExactK_ == true) {
        this->log_.info("calculate K using DEBUG engine.");
        this->getK_exact(P, pK);
    } else {
        this->getK_integralDriven(P, pK);
    }

    //this->cutoffReport();
}


void DfEriX::getK_exact(const TlSymmetricMatrix& P, TlSymmetricMatrix* pK)
{
    assert(pK != NULL);
    // const index_type numOfAOs = this->m_nNumOfAOs;
    pK->resize(this->m_nNumOfAOs);

    const TlOrbitalInfo orbitalInfo((*(this->pPdfParam_))["coordinates"],
                                    (*(this->pPdfParam_))["basis_sets"]);
    const ShellArrayTable shellArrayTable = this->makeShellArrayTable(orbitalInfo);
    const ShellPairArrayTable shellPairArrayTable = this->getShellPairArrayTable(shellArrayTable);
    
    DfEriEngine engine;
    engine.setPrimitiveLevelThreshold(0.0);

    for (int shellTypeP = DfEriX::MAX_SHELL_TYPE -1; shellTypeP >= 0; --shellTypeP) {
        const int maxStepsP = 2 * shellTypeP + 1;
        const ShellArray shellArrayP = shellArrayTable[shellTypeP];
        ShellArray::const_iterator pItEnd = shellArrayP.end();

        for (int shellTypeQ = DfEriX::MAX_SHELL_TYPE -1; shellTypeQ >= 0; --shellTypeQ) {
            const int maxStepsQ = 2 * shellTypeQ + 1;
            const ShellArray shellArrayQ = shellArrayTable[shellTypeQ];
            ShellArray::const_iterator qItEnd = shellArrayQ.end();

            const DfEriEngine::Query queryPQ(0, 0, shellTypeP, shellTypeQ);

            for (int shellTypeR = DfEriX::MAX_SHELL_TYPE -1; shellTypeR >= 0; --shellTypeR) {
                const int maxStepsR = 2 * shellTypeR + 1;
                const ShellArray shellArrayR = shellArrayTable[shellTypeR];
                ShellArray::const_iterator rItEnd = shellArrayR.end();

                for (int shellTypeS = DfEriX::MAX_SHELL_TYPE -1; shellTypeS >= 0; --shellTypeS) {
                    const int maxStepsS = 2 * shellTypeS + 1;
                    const ShellArray shellArrayS = shellArrayTable[shellTypeS];
                    ShellArray::const_iterator sItEnd = shellArrayS.end();

                    const DfEriEngine::Query queryRS(0, 0, shellTypeR, shellTypeS);
        
                    for (ShellArray::const_iterator pIt = shellArrayP.begin(); pIt != pItEnd; ++pIt) {
                        const index_type shellIndexP = *pIt;
                        for (ShellArray::const_iterator qIt = shellArrayQ.begin(); qIt != qItEnd; ++qIt) {
                            const index_type shellIndexQ = *qIt;

                            const DfEriEngine::CGTO_Pair PQ = engine.getCGTO_pair(orbitalInfo,
                                                                                  shellIndexP,
                                                                                  shellIndexQ,
                                                                                  0.0);

                            for (ShellArray::const_iterator rIt = shellArrayR.begin(); rIt != rItEnd; ++rIt) {
                                const index_type shellIndexR = *rIt;
                                for (ShellArray::const_iterator sIt = shellArrayS.begin(); sIt != sItEnd; ++sIt) {
                                    const index_type shellIndexS = *sIt;

                                    const DfEriEngine::CGTO_Pair RS = engine.getCGTO_pair(orbitalInfo,
                                                                                          shellIndexR,
                                                                                          shellIndexS,
                                                                                          0.0);
                                    engine.calc(queryPQ, queryRS, PQ, RS);

                                    int index = 0;
                                    for (int i = 0; i < maxStepsP; ++i) {
                                        const int indexP = shellIndexP + i;
                                        for (int j = 0; j < maxStepsQ; ++j) {
                                            const int indexQ = shellIndexQ + j;
                                            
                                            for (int k = 0; k < maxStepsR; ++k) {
                                                const int indexR = shellIndexR + k;
                                                for (int l = 0; l < maxStepsS; ++l) {
                                                    const int indexS = shellIndexS + l;
                                                    
                                                    if (indexP >= indexR) {
                                                        const double P_qs = P.get(indexQ, indexS);
                                                        const double value = engine.WORK[index];
                                                        pK->add(indexP, indexR, -1.0 * P_qs * value);
                                                    }
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

}


void DfEriX::getK_integralDriven(const TlSymmetricMatrix& P, TlSymmetricMatrix* pK)
{
    assert(pK != NULL);
    TlTime time_all;
    time_all.start();

    const index_type numOfAOs = this->m_nNumOfAOs;
    pK->resize(numOfAOs);
    
    const TlOrbitalInfo orbitalInfo((*(this->pPdfParam_))["coordinates"],
                                    (*(this->pPdfParam_))["basis_sets"]);
    // const ShellArrayTable shellArrayTable = this->makeShellArrayTable(orbitalInfo);

    // ShellPairArrayTable shellPairArrayTable = this->getShellPairArrayTable(shellArrayTable);
    // shellPairArrayTable = this->selectShellPairArrayTableByDensity(shellPairArrayTable,
    //                                                                orbitalInfo);
    
    const TlSparseSymmetricMatrix schwarzTable = this->makeSchwarzTable(orbitalInfo);

#ifdef DEBUG_K
    // for debug
    this->IA_K_ID1_.resize(numOfAOs);
    this->IA_K_ID2_.resize(numOfAOs);
    this->IA_K_ID3_.resize(numOfAOs);
    this->IA_K_ID4_.resize(numOfAOs);
#endif // DEBUG_K

    this->createEngines();
    DfTaskCtrl* pDfTaskCtrl = this->getDfTaskCtrlObject();
    pDfTaskCtrl->setCutoffThreshold(this->cutoffThreshold_);
    pDfTaskCtrl->setCutoffEpsilon_density(this->cutoffEpsilon_density_);
    pDfTaskCtrl->setCutoffEpsilon_distribution(this->cutoffEpsilon_distribution_);

    std::vector<DfTaskCtrl::Task4> taskList;
    bool hasTask = pDfTaskCtrl->getQueue4(orbitalInfo,
                                          schwarzTable,
                                          this->grainSize_,
                                          &taskList,
                                          true);
    while (hasTask == true) {
        this->getK_integralDriven_part(orbitalInfo,
                                       taskList,
                                       P, pK);
        hasTask = pDfTaskCtrl->getQueue4(orbitalInfo,
                                         schwarzTable,
                                         this->grainSize_,
                                         &taskList);
    }

    this->finalize(pK);

    delete pDfTaskCtrl;
    pDfTaskCtrl = NULL;
    const double time_ERI_all = this->destroyEngines();
    
    // debug
#ifdef DEBUG_K
    if (this->isDebugOutK_ == true) {
        this->debugoutK_integralDriven();
    }
#endif // DEBUG_K

    //pK->save("K.mat");

    // statics report
    time_all.stop();
    {
        this->log_.info(TlUtils::format("all time: %16.1f sec.", time_all.getElapseTime()));
        this->log_.info(TlUtils::format("ERI time: %16.1f sec.", time_ERI_all));
    }
}


void DfEriX::getK_integralDriven_part(const TlOrbitalInfoObject& orbitalInfo,
                                      const std::vector<DfTaskCtrl::Task4>& taskList,
                                      const TlMatrixObject& P, TlMatrixObject* pK)
{
    const int taskListSize = taskList.size();
    const double pairwisePGTO_cutoffThreshold = this->cutoffEpsilon3_;

#pragma omp parallel
    {
        int threadID = 0;
#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif // _OPENMP

        this->pEriEngines_[threadID].setPrimitiveLevelThreshold(this->cutoffEpsilon3_);

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
                        
            const DfEriEngine::CGTO_Pair PQ = this->pEriEngines_[threadID].getCGTO_pair(orbitalInfo,
                                                                                        shellIndexP,
                                                                                        shellIndexQ,
                                                                                        pairwisePGTO_cutoffThreshold);
            const DfEriEngine::CGTO_Pair RS = this->pEriEngines_[threadID].getCGTO_pair(orbitalInfo,
                                                                                        shellIndexR,
                                                                                        shellIndexS,
                                                                                        pairwisePGTO_cutoffThreshold);
            const DfEriEngine::Query queryPQ(0, 0, shellTypeP, shellTypeQ);
            const DfEriEngine::Query queryRS(0, 0, shellTypeR, shellTypeS);

            this->pEriEngines_[threadID].calc(queryPQ, queryRS, PQ, RS);
                        
            this->storeK_integralDriven(shellIndexP, maxStepsP,
                                        shellIndexQ, maxStepsQ,
                                        shellIndexR, maxStepsR,
                                        shellIndexS, maxStepsS,
                                        this->pEriEngines_[threadID], P, pK);
        }
    }
}


void DfEriX::storeK_integralDriven(const index_type shellIndexP, const int maxStepsP,
                                   const index_type shellIndexQ, const int maxStepsQ,
                                   const index_type shellIndexR, const int maxStepsR,
                                   const index_type shellIndexS, const int maxStepsS,
                                   const DfEriEngine& engine,
                                   const TlMatrixObject& P,
                                   TlMatrixObject* pK)
{
    assert(pK != NULL);

    int index = 0;
    for (int i = 0; i < maxStepsP; ++i) {
        const index_type indexP = shellIndexP + i;
                                    
        for (int j = 0; j < maxStepsQ; ++j) {
            const index_type indexQ = shellIndexQ + j;
            if (indexP < indexQ) {
                // bypass
                index += (maxStepsR * maxStepsS);
                continue;
            }
                                        
            for (int k = 0; k < maxStepsR; ++k) {
                const index_type indexR = shellIndexR + k;
                if (shellIndexQ == shellIndexS) {
                    if (indexP < indexR) {
                        index += maxStepsS;
                        continue;
                    }
                }
                                            
                for (int l = 0; l < maxStepsS; ++l) {
                    const index_type indexS = shellIndexS + l;
                    if (indexS > ((indexP == indexR) ? indexQ : indexR)) {
                        // bypass
                        ++index;
                        continue;
                    }
                                                
                    const double value = engine.WORK[index];
                                                
                    // Eq.1 : (indexP, indexR) <= (indexQ, indexS)
                    const double coefEq1 = ((indexP == indexR) && (indexQ != indexS)) ? 2.0 : 1.0;
                    pK->add(indexP, indexR, - coefEq1 * P.getLocal(indexQ, indexS) * value);
#ifdef DEBUG_K
                    this->IA_K_ID1_.countUp(indexP, indexR, indexQ, indexS, coefEq1);
#endif // DEBUG_K
                    
                    if (indexR != indexS) { // Eq.1 != Eq.2
                        // Eq.2 : (indexP, indexS) <= (indexQ, indexR)
                        const double coefEq2 = ((indexP == indexS) && (indexQ != indexR)) ? 2.0 : 1.0;
                        pK->add(indexP, indexS, - coefEq2 * P.getLocal(indexQ, indexR) * value);
#ifdef DEBUG_K
                        this->IA_K_ID2_.countUp(indexP, indexS, indexQ, indexR, coefEq2);
#endif // DEBUG_K
                    }
                                                
                    if (indexP != indexQ) { // (Eq.1, Eq.2) != (Eq.3, Eq.4)
                        // Eq.3 : (indexQ, indexR) <= (indexP, indexS)
                        if ((indexP != indexR) || (indexQ != indexS)) { // Eq.2 != Eq.3 : !((indexP, indexS) == (indexQ, indexR))
                            const double coefEq3 = ((indexQ == indexR) && (indexP != indexS)) ? 2.0 : 1.0;
                            pK->add(indexQ, indexR, - coefEq3 * P.getLocal(indexP, indexS) * value);
#ifdef DEBUG_K
                            this->IA_K_ID3_.countUp(indexQ, indexR, indexP, indexS, coefEq3);
#endif // DEBUG_K
                        }
                                                    
                        if (indexR != indexS) { // Eq.3 != Eq.4
                            // Eq.4 : (indexQ, indexS) <= (indexP, indexR)
                            const double coefEq4 = ((indexQ == indexS) && (indexP != indexR)) ? 2.0 : 1.0;
                            pK->add(indexQ, indexS, - coefEq4 * P.getLocal(indexP, indexR) * value);
#ifdef DEBUG_K
                            this->IA_K_ID4_.countUp(indexQ, indexS, indexP, indexR, coefEq4);
                            // std::cerr << TlUtils::format("K_EQ4 <%2d %2d %2d %2d>: (%2d %2d %2d %2d), coef=%d, value=%d",
                            //                              shellIndexP, shellIndexQ, shellIndexR, shellIndexS,
                            //                              indexQ, indexS, indexP, indexR,
                            //                              (int)coefEq4,
                            //                              this->IA_K_ID4_.getCount(indexQ, indexS, indexP, indexR))
                            //           << std::endl;
#endif // DEBUG_K
                        }
                    }
                    ++index;
                }
            }
                                        
        }
    }
}


void DfEriX::debugoutK_integralDriven() const 
{
#ifdef DEBUG_K
    const index_type numOfAOs = this->m_nNumOfAOs;
    for (index_type i = 0; i < numOfAOs; ++i) {
        for (index_type j = 0; j <= i; ++j) {
            std::cerr << TlUtils::format(">>>>K(%2d,%2d)", i, j)
                      << std::endl;
            for (index_type k = 0; k < numOfAOs; ++k) {
                for (index_type l = 0; l <= k; ++l) {
                    const int counter1 = this->IA_K_ID1_.getCount(i, j, k, l);
                    const int counter2 = this->IA_K_ID2_.getCount(i, j, k, l);
                    const int counter3 = this->IA_K_ID3_.getCount(i, j, k, l);
                    const int counter4 = this->IA_K_ID4_.getCount(i, j, k, l);
                    const int counter = counter1 + counter2 + counter3 + counter4;
                    std::string YN = "  ";
                    if (counter != ((k == l) ? 1 : 2)) {
                        YN = "NG";
                    }
                    std::cerr << TlUtils::format("K(%2d,%2d) <= (%2d,%2d) %2d(%2d,%2d,%2d,%2d) %s",
                                                 i, j, k, l, counter,
                                                 counter1, counter2, counter3, counter4,
                                                 YN.c_str())
                              << std::endl;
                }
            }
            std::cerr << std::endl;
        }
    }
#endif // DEBUG_K
}


void DfEriX::getForceK(const TlSymmetricMatrix& P, TlMatrix* pForce)
{
    assert(pForce != NULL);
    //this->clearCutoffStats();

    // カットオフ値の設定
    const double maxDeltaP = P.getMaxAbsoluteElement();
    if (maxDeltaP < 1.0) {
        this->cutoffThreshold_ /= std::fabs(maxDeltaP);
        this->log_.info(TlUtils::format(" new cutoff threshold = % e",
                                        this->cutoffThreshold_));
    }

    pForce->resize(this->m_nNumOfAtoms, 3);
    
    const TlOrbitalInfo orbitalInfo((*(this->pPdfParam_))["coordinates"],
                                    (*(this->pPdfParam_))["basis_sets"]);
    const TlSparseSymmetricMatrix schwarzTable = this->makeSchwarzTable(orbitalInfo);

    this->createEngines();
    DfTaskCtrl* pDfTaskCtrl = this->getDfTaskCtrlObject();

    std::vector<DfTaskCtrl::Task4> taskList;
    pDfTaskCtrl->setCutoffThreshold(this->cutoffThreshold_);
    pDfTaskCtrl->setCutoffEpsilon_density(this->cutoffEpsilon_density_);
    pDfTaskCtrl->setCutoffEpsilon_distribution(this->cutoffEpsilon_distribution_);

    bool hasTask = pDfTaskCtrl->getQueue_Force4(orbitalInfo,
                                                schwarzTable,
                                                this->grainSize_, &taskList, true);
    while (hasTask == true) {
        this->getForceK_part(orbitalInfo, taskList,
                             P, pForce);
        hasTask = pDfTaskCtrl->getQueue_Force4(orbitalInfo,
                                               schwarzTable,
                                               this->grainSize_, &taskList);
    }

    this->finalize(pForce);

    delete pDfTaskCtrl;
    pDfTaskCtrl = NULL;
    this->destroyEngines();
}


void DfEriX::getForceK_part(const TlOrbitalInfoObject& orbitalInfo,
                            const std::vector<DfTaskCtrl::Task4>& taskList,
                            const TlMatrixObject& P, TlMatrixObject* pForce)
{
    const int taskListSize = taskList.size();
    const double pairwisePGTO_cutoffThreshold = this->cutoffEpsilon3_;
    
#pragma omp parallel
    {
        int threadID = 0;
#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif // _OPENMP
        
        this->pEriEngines_[threadID].setPrimitiveLevelThreshold(this->cutoffEpsilon3_);

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
            const index_type atomIndexA = orbitalInfo.getAtomIndex(shellIndexP);
            const index_type atomIndexB = orbitalInfo.getAtomIndex(shellIndexQ);
            const index_type atomIndexC = orbitalInfo.getAtomIndex(shellIndexR);
            const index_type atomIndexD = orbitalInfo.getAtomIndex(shellIndexS);

            // std::cerr << TlUtils::format("[%d %d|%d %d]",
            //                              shellIndexP, shellIndexQ, shellIndexR, shellIndexS)
            //           << std::endl;
            
            if ((atomIndexA == atomIndexB) && (atomIndexB == atomIndexC) &&
                (atomIndexC == atomIndexD) && (atomIndexD == atomIndexA)) {
                continue;
            }
            
            const DfEriEngine::Query queryRS00(0, 0, shellTypeR, shellTypeS);
            const DfEriEngine::Query queryPQ00(0, 0, shellTypeP, shellTypeQ);
            const DfEriEngine::CGTO_Pair PQ =
                this->pEriEngines_[threadID].getCGTO_pair(orbitalInfo,
                                                          shellIndexP,
                                                          shellIndexQ,
                                                          pairwisePGTO_cutoffThreshold);
            const DfEriEngine::CGTO_Pair RS =
                this->pEriEngines_[threadID].getCGTO_pair(orbitalInfo,
                                                          shellIndexR,
                                                          shellIndexS,
                                                          pairwisePGTO_cutoffThreshold);
            
            this->pEriEngines_[threadID].calcGrad(queryPQ00, queryRS00, PQ, RS);
                        
            this->storeForceK_integralDriven(atomIndexA, atomIndexB,
                                             atomIndexC, atomIndexD,
                                             shellIndexP, maxStepsP,
                                             shellIndexQ, maxStepsQ,
                                             shellIndexR, maxStepsR,
                                             shellIndexS, maxStepsS,
                                             this->pEriEngines_[threadID], P, pForce);
        }
    }
}


DfEriX::ShellArrayTable DfEriX::makeShellArrayTable(const TlOrbitalInfoObject& orbitalInfo)
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


void DfEriX::storeForceK_integralDriven(const int atomIndexA, const int atomIndexB,
                                        const int atomIndexC, const int atomIndexD,
                                        const index_type shellIndexP, const int maxStepsP,
                                        const index_type shellIndexQ, const int maxStepsQ,
                                        const index_type shellIndexR, const int maxStepsR,
                                        const index_type shellIndexS, const int maxStepsS,
                                        const DfEriEngine& engine,
                                        const TlMatrixObject& P,
                                        TlMatrixObject* pForce)
{
    int index = 0;
    this->storeForceK_integralDriven(atomIndexA, atomIndexB,
                                     atomIndexC, atomIndexD,
                                     shellIndexP, maxStepsP,
                                     shellIndexQ, maxStepsQ,
                                     shellIndexR, maxStepsR,
                                     shellIndexS, maxStepsS,
                                     engine, P, pForce,
                                     X, &index);
    this->storeForceK_integralDriven(atomIndexA, atomIndexB,
                                     atomIndexC, atomIndexD,
                                     shellIndexP, maxStepsP,
                                     shellIndexQ, maxStepsQ,
                                     shellIndexR, maxStepsR,
                                     shellIndexS, maxStepsS,
                                     engine, P, pForce,
                                     Y, &index);
    this->storeForceK_integralDriven(atomIndexA, atomIndexB,
                                     atomIndexC, atomIndexD,
                                     shellIndexP, maxStepsP,
                                     shellIndexQ, maxStepsQ,
                                     shellIndexR, maxStepsR,
                                     shellIndexS, maxStepsS,
                                     engine, P, pForce,
                                     Z, &index);
}


void DfEriX::storeForceK_integralDriven(const int atomIndexA, const int atomIndexB,
                                        const int atomIndexC, const int atomIndexD,
                                        const index_type shellIndexP, const int maxStepsP,
                                        const index_type shellIndexQ, const int maxStepsQ,
                                        const index_type shellIndexR, const int maxStepsR,
                                        const index_type shellIndexS, const int maxStepsS,
                                        const DfEriEngine& engine,
                                        const TlMatrixObject& P,
                                        TlMatrixObject* pForce,
                                        const int target, int* pIndex)
{
    for (int stepP = 0; stepP < maxStepsP; ++stepP) {
        const index_type indexP = shellIndexP + stepP;
        const index_type iw = indexP * (indexP -1) / 2;

        for (int stepQ = 0; stepQ < maxStepsQ; ++stepQ) {
            const index_type indexQ = shellIndexQ + stepQ;

            if (indexQ <= indexP) {
                const double cij = (indexP != indexQ) ? 2.0 : 1.0;
                const index_type ij = iw + indexQ;
                
                for (int stepR = 0; stepR < maxStepsR; ++stepR) {
                    const index_type indexR = shellIndexR + stepR;
                    const index_type kw = indexR * (indexR -1) / 2;
                    
                    const double dcik = P.get(indexP, indexR);
                    const double dcjk = P.get(indexQ, indexR);
                    
                    for (int stepS = 0; stepS < maxStepsS; ++stepS) {
                        const index_type indexS = shellIndexS + stepS;
                        const index_type kl = kw + indexS;
                        
                        if ((indexS <= indexR) && (ij >= kl)) {
                            double cijkl = cij * ((indexR != indexS) ? 2.0 : 1.0);
                            cijkl *= (ij != kl) ? 2.0 : 1.0;

                            const double coef = 0.5 * cijkl * (  dcik * P.get(indexQ, indexS)
                                                               + dcjk * P.get(indexP, indexS));
                            const double gradIntA = engine.WORK_A[*pIndex];
                            const double gradIntB = engine.WORK_B[*pIndex];
                            const double gradIntC = engine.WORK_C[*pIndex];
                            const double gradIntD = - (gradIntA + gradIntB + gradIntC);
                            
                            pForce->add(atomIndexA, target, coef * gradIntA);
                            pForce->add(atomIndexB, target, coef * gradIntB);
                            pForce->add(atomIndexC, target, coef * gradIntC);
                            pForce->add(atomIndexD, target, coef * gradIntD);
                        }
                        
                        ++(*pIndex);
                    }
                }
            } else {
                *pIndex += (maxStepsR * maxStepsS);
            }
        }
    }
}

DfEriX::ShellPairArrayTable DfEriX::getShellPairArrayTable(const ShellArrayTable& shellArrayTable)
{
    ShellPairArrayTable shellPairArrayTable(DfEriX::MAX_SHELL_TYPE * DfEriX::MAX_SHELL_TYPE);

    for (int shellTypeP = DfEriX::MAX_SHELL_TYPE -1; shellTypeP >= 0; --shellTypeP) {
        const ShellArray& shellArrayP = shellArrayTable[shellTypeP];
        ShellArray::const_iterator pItEnd = shellArrayP.end();

        for (int shellTypeR = DfEriX::MAX_SHELL_TYPE -1; shellTypeR >= 0; --shellTypeR) {
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


// J. Chem. Phys.,105,2726 (1996)
// eq.32
// DfEriX::ShellArray DfEriX::selectShellArrayByDistribution(const ShellArray& inShellArray,
//                                                           const index_type companionShellIndex,
//                                                           const TlOrbitalInfoObject& orbitalInfo)
// {
//     ShellArray answer;
//     answer.reserve(inShellArray.size());
    
//     const TlPosition posB = orbitalInfo.getPosition(companionShellIndex);
//     const int shellTypeB = orbitalInfo.getShellType(companionShellIndex);
//     // orbitalInfoのPGTOリストは指数が小さい順にソートされているため、
//     // 最初(index=0)の指数のみをチェックすれば良い。
//     const double exponentB = orbitalInfo.getExponent(companionShellIndex, 0);
    
//     // check
//     static const double INV_EQ32_COEF = 1.0 / (std::pow(2.0 * TlMath::PI(), 0.25) * TlMath::PI());
//     const double threshold = this->cutoffEpsilon2_ * INV_EQ32_COEF;
//     ShellArray::const_iterator itEnd = inShellArray.end();
//     for (ShellArray::const_iterator it = inShellArray.begin(); it != itEnd; ++it) {
//         const int shellPairType = orbitalInfo.getShellType(*it) * DfEriX::MAX_SHELL_TYPE + shellTypeB;
//         const double distance2 = posB.squareDistanceFrom(orbitalInfo.getPosition(*it));
//         const double exponentA = orbitalInfo.getExponent(*it, 0);

//         const double zetaP = exponentA + exponentB;
//         const double zeta = exponentA * exponentB / zetaP;

//         const double exponent = - zeta * distance2;
//         const double coef = 1.0 / (std::pow(zetaP, 1.25));

//         if (coef * std::exp(exponent) >= threshold) {
//             answer.push_back(*it);

// #pragma omp atomic
//             ++(this->cutoffAlive_E2_[shellPairType]);
//         }

// #pragma omp atomic
//         ++(this->cutoffAll_E2_[shellPairType]);
//     }

//     // swap technique
//     ShellArray(answer).swap(answer);
    
//     return answer;
// }


// J. Chem. Phys.,105,2726 (1996)
// eq.31
// 1/r cutoff
// DfEriX::ShellPairArrayTable DfEriX::selectShellPairArrayTableByDensity(
//     const ShellPairArrayTable& inShellPairArrayTable,
//     const TlOrbitalInfoObject& orbitalInfo)
// {
//     assert(inShellPairArrayTable.size() == (DfEriX::MAX_SHELL_TYPE * DfEriX::MAX_SHELL_TYPE));
    
//     const double cutoffThreshold = this->cutoffEpsilon1_;
//     const double CONTRIBUTE_COEF = 2.0 * std::pow(TlMath::PI(), 2.5);
    
//     TlFmt& FmT = TlFmt::getInstance();
//     static const double coef[6] = {
//         -0.017450254,
//          0.132520568,
//         -0.047915444,
//          0.792267596,
//         -0.583721015,
//          0.697593555
//     };
//     // l = 1.0のときのgamma値
//     static const double gamma1[6] = {
//         0.03,
//         0.081722098,
//         0.222616709,
//         0.606423483,
//         1.651939972,
//         4.5
//     };
//     std::vector<double> gamma(6);
//     {
//         const double inv_lprime = 1.0 / this->lengthScaleParameter_;
//         const double ll = inv_lprime * inv_lprime;
//         for (int i = 0; i < 6; ++i) {
//             gamma[i] = gamma1[i] * ll;
//         }
//     }

//     const int maxShellPairType = DfEriX::MAX_SHELL_TYPE * DfEriX::MAX_SHELL_TYPE;
//     ShellPairArrayTable answer(maxShellPairType);
//     for (int shellPairType = 0; shellPairType < maxShellPairType; ++shellPairType) {
//         const ShellPairArray& shellPairArray = inShellPairArrayTable[shellPairType];
//         const std::size_t shellPairArraySize = shellPairArray.size();
//         ShellPairArray tmp;
//         tmp.reserve(shellPairArraySize);
        
//         for (std::size_t shellPairIndex = 0; shellPairIndex < shellPairArraySize; ++shellPairIndex) {
//             const index_type shellIndexA = shellPairArray[shellPairIndex].shellIndex1;
//             const index_type shellIndexB = shellPairArray[shellPairIndex].shellIndex2;

//             const TlPosition posA = orbitalInfo.getPosition(shellIndexA);
//             const int numOfContractionsA = orbitalInfo.getCgtoContraction(shellIndexA);
//             const TlPosition posB = orbitalInfo.getPosition(shellIndexB);
//             const int numOfContractionsB = orbitalInfo.getCgtoContraction(shellIndexB);
//             const double AB2 = posB.squareDistanceFrom(posA);
        
//             double judge = 0.0;
//             for (int pgtoIndexA = 0; pgtoIndexA < numOfContractionsA; ++pgtoIndexA) {
//                 const double coefA = orbitalInfo.getCoefficient(shellIndexA, pgtoIndexA);
//                 const double zetaA = orbitalInfo.getExponent(shellIndexA, pgtoIndexA);
                
//                 for (int pgtoIndexB = 0; pgtoIndexB < numOfContractionsB; ++pgtoIndexB) {
//                     const double coefB = orbitalInfo.getCoefficient(shellIndexB, pgtoIndexB);
//                     const double zetaB = orbitalInfo.getExponent(shellIndexB, pgtoIndexB);
                    
//                     const double coefAB = coefA * coefB;
//                     const double zetaAB = zetaA * zetaB;
//                     const double zetaA_B = zetaA + zetaB;
//                     for (int i = 0; i < 6; ++i) {
//                         const double zetaAgamma = zetaA * gamma[i];
//                         const double zetaBgamma = zetaB * gamma[i];
//                         const double param = zetaAB + zetaAgamma + zetaBgamma;
//                         const double term1 = CONTRIBUTE_COEF / (std::sqrt(zetaA_B) * param);
//                         const double term2 = std::exp(-zetaAB * gamma[i] * AB2 / param);
//                         const double T = zetaAB * zetaAB * AB2 / (zetaA_B * param);
//                         double term3 = 0.0;
//                         FmT.getFmT(0, T, &term3);
//                         judge += coefAB * coef[i] * term1 * term2 * term3;
//                     }
//                 }
//             }

//             if (std::fabs(judge) > cutoffThreshold) {
//                 tmp.push_back(shellPairArray[shellPairIndex]);

// #pragma omp atomic
//                 ++(this->cutoffAlive_E1_[shellPairType]);
//             }

// #pragma omp atomic
//             ++(this->cutoffAll_E1_[shellPairType]);
//         }

//         // swap technique
//         ShellPairArray(tmp).swap(tmp);

//         answer[shellPairType] = tmp;
//     }
    
//     return answer;
// }


TlSparseSymmetricMatrix DfEriX::makeSchwarzTable(const TlOrbitalInfoObject& orbitalInfo)
{
    this->log_.info("make Schwartz cutoff table: start");
    const index_type maxShellIndex = orbitalInfo.getNumOfOrbitals();
    TlSparseSymmetricMatrix schwarz(maxShellIndex);

    DfEriEngine engine;
    engine.setPrimitiveLevelThreshold(0.0);
    
    for (index_type shellIndexP = 0; shellIndexP < maxShellIndex; ) {
        const int shellTypeP = orbitalInfo.getShellType(shellIndexP);
        const int maxStepsP = 2 * shellTypeP + 1;

        for (index_type shellIndexQ = 0; shellIndexQ < maxShellIndex; ) {
            const int shellTypeQ = orbitalInfo.getShellType(shellIndexQ);
            const int maxStepsQ = 2 * shellTypeQ + 1;
            
            const DfEriEngine::Query queryPQ(0, 0, shellTypeP, shellTypeQ);
            const DfEriEngine::CGTO_Pair PQ = engine.getCGTO_pair(orbitalInfo,
                                                                  shellIndexP,
                                                                  shellIndexQ,
                                                                  0.0);
            engine.calc(queryPQ, queryPQ, PQ, PQ);

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


// bool DfEriX::isAliveBySchwarzCutoff(const index_type shellIndexP,
//                                     const index_type shellIndexQ,
//                                     const index_type shellIndexR,
//                                     const index_type shellIndexS,
//                                     const int shellQuartetType,
//                                     const TlSparseSymmetricMatrix& schwarzTable,
//                                     const double threshold)
// {
//     bool answer = false;

//     const double sqrt_pqpq = schwarzTable.get(shellIndexP, shellIndexQ);
//     const double sqrt_rsrs = schwarzTable.get(shellIndexR, shellIndexS);

//     if ((sqrt_pqpq * sqrt_rsrs) >= threshold) {
//         answer = true;

// #pragma omp atomic
//         ++(this->cutoffAlive_schwarz_[shellQuartetType]);
//     }

// #pragma omp atomic
//     ++(this->cutoffAll_schwarz_[shellQuartetType]);

//     return answer;
// }

