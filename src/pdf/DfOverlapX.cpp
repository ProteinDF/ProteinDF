#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

#include "DfOverlapX.h"
#include "TlOrbitalInfo.h"
#include "TlOrbitalInfo_Density.h"
#include "TlOrbitalInfo_XC.h"

const int DfOverlapX::MAX_SHELL_TYPE = 2 + 1;

DfOverlapX::DfOverlapX(TlSerializeData* pPdfParam) 
    : DfObject(pPdfParam), pEngines_(NULL)
{
}


DfOverlapX::~DfOverlapX()
{
}


void DfOverlapX::createEngines()
{
    assert(this->pEngines_ == NULL);
    this->log_.info(TlUtils::format("create threads: %d", this->numOfThreads_));
    this->pEngines_ = new DfOverlapEngine[this->numOfThreads_];
}

void DfOverlapX::destroyEngines()
{
    this->log_.info("delete threads");
    if (this->pEngines_ != NULL) {
        delete[] this->pEngines_;
        this->pEngines_ = NULL;
    }
}

DfTaskCtrl* DfOverlapX::getDfTaskCtrlObject() const
{
    DfTaskCtrl* pDfTaskCtrl = new DfTaskCtrl(this->pPdfParam_);
    return pDfTaskCtrl;
}


void DfOverlapX::finalize(TlMatrix* pMtx)
{
    // do nothing
}


void DfOverlapX::finalize(TlSymmetricMatrix* pMtx)
{
    // do nothing
}


void DfOverlapX::finalize(TlVector* pVct)
{
    // do nothing
}


void DfOverlapX::getSpq(TlSymmetricMatrix* pSpq)
{
    assert(pSpq != NULL);
    const index_type numOfAOs = this->m_nNumOfAOs;
    pSpq->resize(numOfAOs);

    const TlOrbitalInfo orbitalInfo((*(this->pPdfParam_))["coordinates"],
                                    (*(this->pPdfParam_))["basis_sets"]);
    this->calcOverlap(orbitalInfo, pSpq);
    this->finalize(pSpq);
}


void DfOverlapX::getSab(TlSymmetricMatrix* pSab)
{
    assert(pSab != NULL);
    const index_type numOfAuxDens = this->m_nNumOfAux;
    pSab->resize(numOfAuxDens);

    const TlOrbitalInfo_Density orbitalInfo_Density((*(this->pPdfParam_))["coordinates"],
                                                    (*(this->pPdfParam_))["basis_sets_j"]);
    this->calcOverlap(orbitalInfo_Density, pSab);
    this->finalize(pSab);
}

void DfOverlapX::getSgd(TlSymmetricMatrix* pSgd)
{
    assert(pSgd != NULL);
    const index_type numOfAuxXC = this->numOfAuxXC_;
    pSgd->resize(numOfAuxXC);

    const TlOrbitalInfo_XC orbitalInfo_XC((*(this->pPdfParam_))["coordinates"],
                                          (*(this->pPdfParam_))["basis_sets_k"]);
    this->calcOverlap(orbitalInfo_XC, pSgd);
    this->finalize(pSgd);
}

void DfOverlapX::getNalpha(TlVector* pNalpha)
{
    assert(pNalpha != NULL);
    const index_type numOfAuxDens = this->m_nNumOfAux;
    pNalpha->resize(numOfAuxDens);
    
    const TlOrbitalInfo_Density orbitalInfo_Density((*(this->pPdfParam_))["coordinates"],
                                                    (*(this->pPdfParam_))["basis_sets_j"]);
    this->calcOverlap(orbitalInfo_Density, pNalpha);
    this->finalize(pNalpha);
}

void DfOverlapX::getOvpMat(const TlOrbitalInfoObject& orbitalInfo,
                           TlSymmetricMatrix* pS)
{
    assert(pS != NULL);
    pS->resize(orbitalInfo.getNumOfOrbitals());

    this->calcOverlap(orbitalInfo, pS);
    this->finalize(pS);
}

void DfOverlapX::getTransMat(const TlOrbitalInfoObject& orbitalInfo1,
                             const TlOrbitalInfoObject& orbitalInfo2,
                             TlMatrix* pTransMat)
{
    assert(pTransMat != NULL);
    pTransMat->resize(orbitalInfo1.getNumOfOrbitals(),
                      orbitalInfo2.getNumOfOrbitals());

    this->calcOverlap(orbitalInfo1, orbitalInfo2, pTransMat);
    this->finalize(pTransMat);
}

void DfOverlapX::get_pqg(const TlVector& myu, TlSymmetricMatrix* pF)
{
    assert(pF != NULL);
    const TlOrbitalInfo orbitalInfo((*(this->pPdfParam_))["coordinates"],
                                    (*(this->pPdfParam_))["basis_sets"]);
    const TlOrbitalInfo_XC orbitalInfo_XC((*(this->pPdfParam_))["coordinates"],
                                          (*(this->pPdfParam_))["basis_sets_k"]);
    pF->resize(orbitalInfo.getNumOfOrbitals());
    this->calcOverlap(orbitalInfo_XC, myu,
                      orbitalInfo, pF);
    this->finalize(pF);
}

void DfOverlapX::get_pqg(const TlVector& myu, const TlVector& eps,
                         TlSymmetricMatrix* pF,
                         TlSymmetricMatrix* pE)
{
    assert(pF != NULL);
    assert(pE != NULL);
    const TlOrbitalInfo orbitalInfo((*(this->pPdfParam_))["coordinates"],
                                    (*(this->pPdfParam_))["basis_sets"]);
    const TlOrbitalInfo_XC orbitalInfo_XC((*(this->pPdfParam_))["coordinates"],
                                          (*(this->pPdfParam_))["basis_sets_k"]);
    pF->resize(orbitalInfo.getNumOfOrbitals());
    this->calcOverlap(orbitalInfo_XC, myu, eps,
                      orbitalInfo, pF, pE);
    this->finalize(pF);
    this->finalize(pE);
}

void DfOverlapX::calcOverlap(const TlOrbitalInfoObject& orbitalInfo,
                             TlMatrixObject* pMatrix)
{
    this->createEngines();
    DfTaskCtrl* pTaskCtrl = this->getDfTaskCtrlObject();

    std::vector<DfTaskCtrl::Task2> taskList;
    bool hasTask = pTaskCtrl->getQueue2(orbitalInfo,
                                       true,
                                       this->grainSize_, &taskList, true);
    while (hasTask == true) {
        this->calcOverlap_part(orbitalInfo,
                               taskList,
                               pMatrix);
        
        hasTask = pTaskCtrl->getQueue2(orbitalInfo,
                                       true,
                                       this->grainSize_, &taskList);
    }

    pTaskCtrl->cutoffReport();
    delete pTaskCtrl;
    pTaskCtrl = NULL;
    this->destroyEngines();
}


void DfOverlapX::calcOverlap(const TlOrbitalInfoObject& orbitalInfo1,
                             const TlOrbitalInfoObject& orbitalInfo2,
                             TlMatrixObject* pMatrix)
{
    this->createEngines();
    DfTaskCtrl* pTaskCtrl = this->getDfTaskCtrlObject();

    std::vector<DfTaskCtrl::Task2> taskList;
    bool hasTask = pTaskCtrl->getQueue2(orbitalInfo1, orbitalInfo2,
                                       true,
                                       this->grainSize_, &taskList, true);
    while (hasTask == true) {
        this->calcOverlap_part(orbitalInfo1,
                               orbitalInfo2,
                               taskList,
                               pMatrix);
        
        hasTask = pTaskCtrl->getQueue2(orbitalInfo1, orbitalInfo2,
                                       true,
                                       this->grainSize_, &taskList);
    }

    pTaskCtrl->cutoffReport();
    delete pTaskCtrl;
    pTaskCtrl = NULL;
    this->destroyEngines();
}


void DfOverlapX::calcOverlap(const TlOrbitalInfoObject& orbitalInfo,
                             TlVectorObject* pVector)
{
    this->createEngines();
    DfTaskCtrl* pTaskCtrl = this->getDfTaskCtrlObject();

    std::vector<DfTaskCtrl::Task> taskList;
    bool hasTask = pTaskCtrl->getQueue(orbitalInfo,
                                       this->grainSize_, &taskList, true);
    while (hasTask == true) {
        this->calcOverlap_part(orbitalInfo,
                               taskList,
                               pVector);
        
        hasTask = pTaskCtrl->getQueue(orbitalInfo,
                                      this->grainSize_, &taskList);
    }

    pTaskCtrl->cutoffReport();
    delete pTaskCtrl;
    pTaskCtrl = NULL;
    this->destroyEngines();
}

void DfOverlapX::calcOverlap(const TlOrbitalInfoObject& orbitalInfo_XC,
                             const TlVector& myu,
                             const TlOrbitalInfoObject& orbitalInfo,
                             TlMatrixObject* pF)
{
    this->createEngines();
    DfTaskCtrl* pTaskCtrl = this->getDfTaskCtrlObject();

    std::vector<DfTaskCtrl::Task2> taskList;
    bool hasTask = pTaskCtrl->getQueue2(orbitalInfo,
                                       true,
                                       this->grainSize_, &taskList, true);
    while (hasTask == true) {
        this->calcOverlap_part(orbitalInfo_XC,
                               myu,
                               orbitalInfo,
                               taskList,
                               pF);
        
        hasTask = pTaskCtrl->getQueue2(orbitalInfo,
                                       true,
                                       this->grainSize_, &taskList);
    }

    pTaskCtrl->cutoffReport();
    delete pTaskCtrl;
    pTaskCtrl = NULL;
    this->destroyEngines();
}

void DfOverlapX::calcOverlap(const TlOrbitalInfoObject& orbitalInfo_XC,
                             const TlVector& myu,
                             const TlVector& eps,
                             const TlOrbitalInfoObject& orbitalInfo,
                             TlMatrixObject* pF,
                             TlMatrixObject* pE)
{
    this->createEngines();
    DfTaskCtrl* pTaskCtrl = this->getDfTaskCtrlObject();

    std::vector<DfTaskCtrl::Task2> taskList;
    bool hasTask = pTaskCtrl->getQueue2(orbitalInfo,
                                       true,
                                       this->grainSize_, &taskList, true);
    while (hasTask == true) {
        this->calcOverlap_part(orbitalInfo_XC,
                               myu, eps,
                               orbitalInfo,
                               taskList,
                               pF, pE);
        
        hasTask = pTaskCtrl->getQueue2(orbitalInfo,
                                       true,
                                       this->grainSize_, &taskList);
    }

    pTaskCtrl->cutoffReport();
    delete pTaskCtrl;
    pTaskCtrl = NULL;
    this->destroyEngines();
}

void DfOverlapX::calcOverlap_part(const TlOrbitalInfoObject& orbitalInfo,
                                  const std::vector<DfTaskCtrl::Task2>& taskList,
                                  TlMatrixObject* pMatrix)
{
    // 第三、四中心点の固定値を設定
    // static const TlPosition posR(0.0, 0.0, 0.0);
    // static DfOverlapEngine::PGTO pgtoR(1.0, 0.0);
    // static DfOverlapEngine::PGTOs pgtosR(1);
    // pgtosR[0] = pgtoR;
    // static const TlPosition posS(0.0, 0.0, 0.0);
    // static DfOverlapEngine::PGTO pgtoS(1.0, 0.0);
    // static DfOverlapEngine::PGTOs pgtosS(1);
    // pgtosS[0] = pgtoS;

    const int taskListSize = taskList.size();

#pragma omp parallel
    {
        int threadID = 0;
#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif // _OPENMP
        assert(threadID < this->numOfThreads_);

#pragma omp for schedule(runtime)
        for (int i = 0; i < taskListSize; ++i) {
            const index_type shellIndexP = taskList[i].shellIndex1;
            const index_type shellIndexQ = taskList[i].shellIndex2;
            
            const int shellTypeP = orbitalInfo.getShellType(shellIndexP);
            const int shellTypeQ = orbitalInfo.getShellType(shellIndexQ);
            const int maxStepsP = 2 * shellTypeP + 1;
            const int maxStepsQ = 2 * shellTypeQ + 1;
            
            this->pEngines_[threadID].calc(0, orbitalInfo, shellIndexP,
                                           0, orbitalInfo, shellIndexQ,
                                           0, orbitalInfo, -1,
                                           0, orbitalInfo, -1);
            
            int index = 0;
            for (int stepP = 0; stepP < maxStepsP; ++stepP) {
                const index_type globalShellIndexP = shellIndexP + stepP;

                for (int stepQ = 0; stepQ < maxStepsQ; ++stepQ) {
                    const index_type globalShellIndexQ = shellIndexQ + stepQ;

                    if ((shellIndexP != shellIndexQ) || (globalShellIndexP >= globalShellIndexQ)) {
                        pMatrix->add(globalShellIndexP, globalShellIndexQ, this->pEngines_[threadID].WORK[index]);
                    }
                    ++index;
                }
            }
        }
    }
}

void DfOverlapX::calcOverlap_part(const TlOrbitalInfoObject& orbitalInfo1,
                                  const TlOrbitalInfoObject& orbitalInfo2,
                                  const std::vector<DfTaskCtrl::Task2>& taskList,
                                  TlMatrixObject* pMatrix)
{
    // 第三、四中心点の固定値を設定
    // static const TlPosition posR(0.0, 0.0, 0.0);
    // static DfOverlapEngine::PGTO pgtoR(1.0, 0.0);
    // static DfOverlapEngine::PGTOs pgtosR(1);
    // pgtosR[0] = pgtoR;
    // static const TlPosition posS(0.0, 0.0, 0.0);
    // static DfOverlapEngine::PGTO pgtoS(1.0, 0.0);
    // static DfOverlapEngine::PGTOs pgtosS(1);
    // pgtosS[0] = pgtoS;

    const int taskListSize = taskList.size();

#pragma omp parallel
    {
        int threadID = 0;
#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif // _OPENMP
        assert(threadID < this->numOfThreads_);

#pragma omp for schedule(runtime)
        for (int i = 0; i < taskListSize; ++i) {
            const index_type shellIndexP = taskList[i].shellIndex1;
            const index_type shellIndexQ = taskList[i].shellIndex2;
            
            const int shellTypeP = orbitalInfo1.getShellType(shellIndexP);
            const int shellTypeQ = orbitalInfo2.getShellType(shellIndexQ);
            const int maxStepsP = 2 * shellTypeP + 1;
            const int maxStepsQ = 2 * shellTypeQ + 1;
            // const TlPosition posP = orbitalInfo1.getPosition(shellIndexP);
            // const TlPosition posQ = orbitalInfo2.getPosition(shellIndexQ);
            // const DfOverlapEngine::PGTOs pgtosP = DfOverlapEngine::getPGTOs(orbitalInfo1, shellIndexP);
            // const DfOverlapEngine::PGTOs pgtosQ = DfOverlapEngine::getPGTOs(orbitalInfo2, shellIndexQ);
            // const DfOverlapEngine::Query query(0, 0, 0, 0, shellTypeP, shellTypeQ, 0, 0);
            
            // this->pEngines_[threadID].calc0(query, posP, posQ, posR, posS, pgtosP, pgtosQ, pgtosR, pgtosS);
            this->pEngines_[threadID].calc(0, orbitalInfo1, shellIndexP,
                                           0, orbitalInfo2, shellIndexQ,
                                           0, orbitalInfo1, -1,
                                           0, orbitalInfo1, -1);
            
            int index = 0;
            for (int stepP = 0; stepP < maxStepsP; ++stepP) {
                const index_type globalShellIndexP = shellIndexP + stepP;

                for (int stepQ = 0; stepQ < maxStepsQ; ++stepQ) {
                    const index_type globalShellIndexQ = shellIndexQ + stepQ;

                    pMatrix->add(globalShellIndexP, globalShellIndexQ, this->pEngines_[threadID].WORK[index]);
                    ++index;
                }
            }
        }
    }
}

void DfOverlapX::calcOverlap_part(const TlOrbitalInfoObject& orbitalInfo_XC,
                                  const TlVector& myu,
                                  const TlOrbitalInfoObject& orbitalInfo,
                                  const std::vector<DfTaskCtrl::Task2>& taskList,
                                  TlMatrixObject* pMatrix)
{
    // 第四中心点の固定値を設定
    // static const TlPosition posS(0.0, 0.0, 0.0);
    // static DfOverlapEngine::PGTO pgtoS(1.0, 0.0);
    // static DfOverlapEngine::PGTOs pgtosS(1);
    // pgtosS[0] = pgtoS;

    const index_type numOfAuxXC = orbitalInfo_XC.getNumOfOrbitals();
    const int taskListSize = taskList.size();

#pragma omp parallel
    {
        int threadID = 0;
#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif // _OPENMP

#pragma omp for schedule(runtime)
        for (int i = 0; i < taskListSize; ++i) {
            const index_type shellIndexP = taskList[i].shellIndex1;
            const index_type shellIndexQ = taskList[i].shellIndex2;
            
            const int shellTypeP = orbitalInfo.getShellType(shellIndexP);
            const int shellTypeQ = orbitalInfo.getShellType(shellIndexQ);
            const int maxStepsP = 2 * shellTypeP + 1;
            const int maxStepsQ = 2 * shellTypeQ + 1;
            // const TlPosition posP = orbitalInfo.getPosition(shellIndexP);
            // const TlPosition posQ = orbitalInfo.getPosition(shellIndexQ);
            // const DfOverlapEngine::PGTOs pgtosP = DfOverlapEngine::getPGTOs(orbitalInfo, shellIndexP);
            // const DfOverlapEngine::PGTOs pgtosQ = DfOverlapEngine::getPGTOs(orbitalInfo, shellIndexQ);

            for (index_type shellIndexR = 0; shellIndexR < numOfAuxXC; ) {
                const int shellTypeR = orbitalInfo_XC.getShellType(shellIndexR);
                const int maxStepsR = 2 * shellTypeR + 1;
                // const TlPosition posR = orbitalInfo_XC.getPosition(shellIndexR);
                // const DfOverlapEngine::PGTOs pgtosR = DfOverlapEngine::getPGTOs(orbitalInfo_XC, shellIndexR);

                // const DfOverlapEngine::Query query(0, 0, 0, 0, shellTypeP, shellTypeQ, shellTypeR, 0);
            
                // this->pEngines_[threadID].calc0(query, posP, posQ, posR, posS, pgtosP, pgtosQ, pgtosR, pgtosS);
                this->pEngines_[threadID].calc(0, orbitalInfo, shellIndexP,
                                               0, orbitalInfo, shellIndexQ,
                                               0, orbitalInfo_XC, shellIndexR,
                                               0, orbitalInfo_XC, -1);
                
                int index = 0;
                for (int stepP = 0; stepP < maxStepsP; ++stepP) {
                    const index_type globalShellIndexP = shellIndexP + stepP;
                    
                    for (int stepQ = 0; stepQ < maxStepsQ; ++stepQ) {
                        const index_type globalShellIndexQ = shellIndexQ + stepQ;
                        
                        for (int stepR = 0; stepR < maxStepsR; ++stepR) {
                            const index_type globalShellIndexR = shellIndexR + stepR;
                            if ((shellIndexP != shellIndexQ) || (globalShellIndexP >= globalShellIndexQ)) {
                                pMatrix->add(globalShellIndexP, globalShellIndexQ, myu.get(globalShellIndexR) * this->pEngines_[threadID].WORK[index]);
                            }
                        }
                        ++index;
                    }
                }
                
                shellIndexR += maxStepsR;
            }
        }
    }
}

void DfOverlapX::calcOverlap_part(const TlOrbitalInfoObject& orbitalInfo_XC,
                                  const TlVector& myu,
                                  const TlVector& eps,
                                  const TlOrbitalInfoObject& orbitalInfo,
                                  const std::vector<DfTaskCtrl::Task2>& taskList,
                                  TlMatrixObject* pF,
                                  TlMatrixObject* pE)
{
    // 第四中心点の固定値を設定
    // static const TlPosition posS(0.0, 0.0, 0.0);
    // static DfOverlapEngine::PGTO pgtoS(1.0, 0.0);
    // static DfOverlapEngine::PGTOs pgtosS(1);
    // pgtosS[0] = pgtoS;

    const index_type numOfAuxXC = orbitalInfo_XC.getNumOfOrbitals();
    const int taskListSize = taskList.size();

#pragma omp parallel
    {
        int threadID = 0;
#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif // _OPENMP

#pragma omp for schedule(runtime)
        for (int i = 0; i < taskListSize; ++i) {
            const index_type shellIndexP = taskList[i].shellIndex1;
            const index_type shellIndexQ = taskList[i].shellIndex2;
            
            const int shellTypeP = orbitalInfo.getShellType(shellIndexP);
            const int shellTypeQ = orbitalInfo.getShellType(shellIndexQ);
            const int maxStepsP = 2 * shellTypeP + 1;
            const int maxStepsQ = 2 * shellTypeQ + 1;
            // const TlPosition posP = orbitalInfo.getPosition(shellIndexP);
            // const TlPosition posQ = orbitalInfo.getPosition(shellIndexQ);
            // const DfOverlapEngine::PGTOs pgtosP = DfOverlapEngine::getPGTOs(orbitalInfo, shellIndexP);
            // const DfOverlapEngine::PGTOs pgtosQ = DfOverlapEngine::getPGTOs(orbitalInfo, shellIndexQ);

            for (index_type shellIndexR = 0; shellIndexR < numOfAuxXC; ) {
                const int shellTypeR = orbitalInfo_XC.getShellType(shellIndexR);
                const int maxStepsR = 2 * shellTypeR + 1;
                // const TlPosition posR = orbitalInfo_XC.getPosition(shellIndexR);
                // const DfOverlapEngine::PGTOs pgtosR = DfOverlapEngine::getPGTOs(orbitalInfo_XC, shellIndexR);

                // const DfOverlapEngine::Query query(0, 0, 0, 0, shellTypeP, shellTypeQ, shellTypeR, 0);
            
                // this->pEngines_[threadID].calc0(query, posP, posQ, posR, posS, pgtosP, pgtosQ, pgtosR, pgtosS);
                this->pEngines_[threadID].calc(0, orbitalInfo, shellIndexP,
                                               0, orbitalInfo, shellIndexQ,
                                               0, orbitalInfo_XC, shellIndexR,
                                               0, orbitalInfo_XC, -1);
                
                int index = 0;
                for (int stepP = 0; stepP < maxStepsP; ++stepP) {
                    const index_type globalShellIndexP = shellIndexP + stepP;
                    
                    for (int stepQ = 0; stepQ < maxStepsQ; ++stepQ) {
                        const index_type globalShellIndexQ = shellIndexQ + stepQ;
                        
                        for (int stepR = 0; stepR < maxStepsR; ++stepR) {
                            const index_type globalShellIndexR = shellIndexR + stepR;
                            if ((shellIndexP != shellIndexQ) || (globalShellIndexP >= globalShellIndexQ)) {
                                const double value = this->pEngines_[threadID].WORK[index];
                                pF->add(globalShellIndexP, globalShellIndexQ, myu.get(globalShellIndexR) * value);
                                pE->add(globalShellIndexP, globalShellIndexQ, eps.get(globalShellIndexR) * value);
                            }
                        }
                        ++index;
                    }
                }
                
                shellIndexR += maxStepsR;
            }
        }
    }
}

void DfOverlapX::calcOverlap_part(const TlOrbitalInfoObject& orbitalInfo,
                                  const std::vector<DfTaskCtrl::Task>& taskList,
                                  TlVectorObject* pVector)
{
    // 第2, 第3, 第4中心点の固定値を設定
    // static const TlPosition posQ(0.0, 0.0, 0.0);
    // static DfOverlapEngine::PGTO pgtoQ(1.0, 0.0);
    // static DfOverlapEngine::PGTOs pgtosQ(1);
    // pgtosQ[0] = pgtoQ;
    // static const TlPosition posR(0.0, 0.0, 0.0);
    // static DfOverlapEngine::PGTO pgtoR(1.0, 0.0);
    // static DfOverlapEngine::PGTOs pgtosR(1);
    // pgtosR[0] = pgtoR;
    // static const TlPosition posS(0.0, 0.0, 0.0);
    // static DfOverlapEngine::PGTO pgtoS(1.0, 0.0);
    // static DfOverlapEngine::PGTOs pgtosS(1);
    // pgtosS[0] = pgtoS;

    const int taskListSize = taskList.size();

#pragma omp parallel
    {
        int threadID = 0;
#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif // _OPENMP

#pragma omp for schedule(runtime)
        for (int i = 0; i < taskListSize; ++i) {
            const index_type shellIndexP = taskList[i].shellIndex1;
            
            const int shellTypeP = orbitalInfo.getShellType(shellIndexP);
            const int maxStepsP = 2 * shellTypeP + 1;
            // const TlPosition posP = orbitalInfo.getPosition(shellIndexP);
            // const DfOverlapEngine::PGTOs pgtosP = DfOverlapEngine::getPGTOs(orbitalInfo, shellIndexP);
            // const DfOverlapEngine::Query query(0, 0, 0, 0, shellTypeP, 0, 0, 0);
            
            // this->pEngines_[threadID].calc0(query, posP, posQ, posR, posS, pgtosP, pgtosQ, pgtosR, pgtosS);
            this->pEngines_[threadID].calc(0, orbitalInfo, shellIndexP,
                                           0, orbitalInfo, -1,
                                           0, orbitalInfo, -1,
                                           0, orbitalInfo, -1);
            
            int index = 0;
            for (int stepP = 0; stepP < maxStepsP; ++stepP) {
                const index_type globalShellIndexP = shellIndexP + stepP;

                pVector->add(globalShellIndexP, this->pEngines_[threadID].WORK[index]);
                ++index;
            }
        }
    }
}


void DfOverlapX::getForce(const TlSymmetricMatrix& W,
                          TlMatrix* pForce)
{
    assert(pForce != NULL);
    pForce->resize(this->m_nNumOfAtoms, 3);
    
    const TlOrbitalInfo orbitalInfo((*(this->pPdfParam_))["coordinates"],
                                    (*(this->pPdfParam_))["basis_sets"]);
    const ShellArrayTable shellArrayTable = this->makeShellArrayTable(orbitalInfo);

    this->createEngines();
    
    for (int shellTypeP = DfOverlapX::MAX_SHELL_TYPE -1; shellTypeP >= 0; --shellTypeP) {
        const ShellArray shellArrayP = shellArrayTable[shellTypeP];
        const index_type shellArraySizeP = shellArrayP.size();
        
        for (int shellTypeQ = DfOverlapX::MAX_SHELL_TYPE -1; shellTypeQ >= 0; --shellTypeQ) {
            const ShellArray shellArrayQ = shellArrayTable[shellTypeQ];
            // const index_type shellArraySizeQ = shellArrayQ.size();
            
            for (index_type p = 0; p < shellArraySizeP; ++p) {
                const index_type shellIndexP = shellArrayP[p];

                this->getForce_partProc(orbitalInfo,
                                        shellTypeP, shellTypeQ,
                                        shellIndexP,
                                        shellArrayQ,
                                        W, pForce);
            }
        }
    }

    this->destroyEngines();
}


void DfOverlapX::getForce_partProc(const TlOrbitalInfoObject& orbitalInfo,
                                   const int shellTypeP, const int shellTypeQ,
                                   const index_type shellIndexP,
                                   const ShellArray& shellArrayQ,
                                   const TlSymmetricMatrix& W,
                                   TlMatrix* pForce)
{
    // 第3,第4中心点の固定値を設定
    // static const TlPosition posR(0.0, 0.0, 0.0);
    // static const DfOverlapEngine::PGTO pgtoR(1.0, 0.0);
    // static DfOverlapEngine::PGTOs pgtosR(1);
    // pgtosR[0] = pgtoR;
    // static const TlPosition posS(0.0, 0.0, 0.0);
    // static const DfOverlapEngine::PGTO pgtoS(1.0, 0.0);
    // static DfOverlapEngine::PGTOs pgtosS(1);
    // pgtosS[0] = pgtoS;

    const int maxStepsP = 2 * shellTypeP + 1;
    const int maxStepsQ = 2 * shellTypeQ + 1;
    const std::size_t shellArraySizeQ = shellArrayQ.size();
    // const DfOverlapEngine::Query query(1, 0, 0, 0, shellTypeP, shellTypeQ, 0, 0);
    
    // const DfOverlapEngine::PGTOs pgtosP = DfOverlapEngine::getPGTOs(orbitalInfo,
    //                                                                 shellIndexP);
    // const TlPosition posP = orbitalInfo.getPosition(shellIndexP);
    const int atomIndexA = orbitalInfo.getAtomIndex(shellIndexP);

#pragma omp parallel
    {
        int threadID = 0;
#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif // _OPENMP

#pragma omp for schedule(runtime)
        for (std::size_t q = 0; q < shellArraySizeQ; ++q) {
            const index_type shellIndexQ = shellArrayQ[q];
            // 重複ペアの排除
//                     if (shellIndexP < shellIndexQ) {
//                         continue;
//                     }
                    
            // const TlPosition posQ = orbitalInfo.getPosition(shellIndexQ);
            // const DfOverlapEngine::PGTOs pgtosQ = DfOverlapEngine::getPGTOs(orbitalInfo,
            //                                                                 shellIndexQ);
            const int atomIndexB = orbitalInfo.getAtomIndex(shellIndexQ);
                
            // this->pEngines_[threadID].calc0(query, posP, posQ, posR, posS, pgtosP, pgtosQ, pgtosR, pgtosS);
            this->pEngines_[threadID].calc(1, orbitalInfo, shellIndexP,
                                           0, orbitalInfo, shellIndexQ,
                                           0, orbitalInfo, -1,
                                           0, orbitalInfo, -1);

            int index = 0;
            for (int stepP = 0; stepP < maxStepsP; ++stepP) {
                const index_type indexP = shellIndexP + stepP;
                for (int stepQ = 0; stepQ < maxStepsQ; ++stepQ) {
                    const index_type indexQ = shellIndexQ + stepQ;

                    double coef = W.get(indexP, indexQ);
                    const double dSdA = this->pEngines_[threadID].WORK[index];
                    const double dSdB = - dSdA;
#pragma omp critical(DfOverlapX__getForce)
                    {
                        pForce->add(atomIndexA, X, coef * dSdA);
                        pForce->add(atomIndexB, X, coef * dSdB);
                    }
                    ++index;
                }
            }
            for (int stepP = 0; stepP < maxStepsP; ++stepP) {
                const index_type indexP = shellIndexP + stepP;
                for (int stepQ = 0; stepQ < maxStepsQ; ++stepQ) {
                    const index_type indexQ = shellIndexQ + stepQ;

                    double coef = W.get(indexP, indexQ);
                    const double dSdA = this->pEngines_[threadID].WORK[index];
                    const double dSdB = - dSdA;
#pragma omp critical(DfOverlapX__getForce)
                    {
                        pForce->add(atomIndexA, Y, coef * dSdA);
                        pForce->add(atomIndexB, Y, coef * dSdB);
                    }
                    ++index;
                }
            }
            for (int stepP = 0; stepP < maxStepsP; ++stepP) {
                const index_type indexP = shellIndexP + stepP;
                for (int stepQ = 0; stepQ < maxStepsQ; ++stepQ) {
                    const index_type indexQ = shellIndexQ + stepQ;

                    double coef = W.get(indexP, indexQ);
                    const double dSdA = this->pEngines_[threadID].WORK[index];
                    const double dSdB = - dSdA;
#pragma omp critical(DfOverlapX__getForce)
                    {
                        pForce->add(atomIndexA, Z, coef * dSdA);
                        pForce->add(atomIndexB, Z, coef * dSdB);
                    }
                    ++index;
                }
            }
        }
    }
}


DfOverlapX::ShellArrayTable DfOverlapX::makeShellArrayTable(const TlOrbitalInfoObject& orbitalInfo)
{
    ShellArrayTable shellArrayTable(MAX_SHELL_TYPE);
    const int maxShellIndex = orbitalInfo.getNumOfOrbitals();

    int shellIndex = 0;
    while (shellIndex < maxShellIndex) {
        // shellType: 0=s, 1=p, 2=d
        const int shellType = orbitalInfo.getShellType(shellIndex);
        const int steps = 2 * shellType +1;

        shellArrayTable[shellType].push_back(shellIndex);
        
        shellIndex += steps;
    }

    return shellArrayTable;
}


// DfOverlapEngine::PGTOs DfOverlapX::getPGTOs(const TlOrbitalInfoObject& orbitalInfo,
//                                             const int shellIndex)
// {
//     const int numOfContractions = orbitalInfo.getCgtoContraction(shellIndex);
//     DfOverlapEngine::PGTOs pgtos(numOfContractions);

//     for (int i = 0; i < numOfContractions; ++i) {
//         const DfOverlapEngine::PGTO pgto(orbitalInfo.getCoefficient(shellIndex, i),
//                                          orbitalInfo.getExponent(shellIndex, i));
//         pgtos[i] = pgto;
//     }
    
//     return pgtos;
// }


void DfOverlapX::getGradient(const TlOrbitalInfoObject& orbitalInfo,
                             TlMatrix* pMatX, TlMatrix* pMatY, TlMatrix* pMatZ)
{
    assert(pMatX != NULL);
    assert(pMatY != NULL);
    assert(pMatZ != NULL);
    
    const index_type numOfAOs = orbitalInfo.getNumOfOrbitals();
    pMatX->resize(numOfAOs, numOfAOs);
    pMatY->resize(numOfAOs, numOfAOs);
    pMatZ->resize(numOfAOs, numOfAOs);

    const ShellArrayTable shellArrayTable = this->makeShellArrayTable(orbitalInfo);

    this->createEngines();
    DfTaskCtrl* pTaskCtrl = this->getDfTaskCtrlObject();

    std::vector<DfTaskCtrl::Task2> taskList;
    bool hasTask = pTaskCtrl->getQueue2(orbitalInfo,
                                       true,
                                       this->grainSize_, &taskList, true);
    while (hasTask == true) {
        this->getGradient_partProc(orbitalInfo,
                                   taskList,
                                   pMatX, pMatY, pMatZ);
        
        hasTask = pTaskCtrl->getQueue2(orbitalInfo,
                                       true,
                                       this->grainSize_, &taskList);
    }

    pTaskCtrl->cutoffReport();
    delete pTaskCtrl;
    pTaskCtrl = NULL;
    this->destroyEngines();

    this->finalize(pMatX);
    this->finalize(pMatY);
    this->finalize(pMatZ);
}


void DfOverlapX::getGradient_partProc(const TlOrbitalInfoObject& orbitalInfo,
                                      const std::vector<DfTaskCtrl::Task2>& taskList,
                                      TlMatrixObject* pMatX, TlMatrixObject* pMatY, TlMatrixObject* pMatZ)
{
    // 第3,第4中心点の固定値を設定
    // static const TlPosition posR(0.0, 0.0, 0.0);
    // static const DfOverlapEngine::PGTO pgtoR(1.0, 0.0);
    // static DfOverlapEngine::PGTOs pgtosR(1);
    // pgtosR[0] = pgtoR;
    // static const TlPosition posS(0.0, 0.0, 0.0);
    // static const DfOverlapEngine::PGTO pgtoS(1.0, 0.0);
    // static DfOverlapEngine::PGTOs pgtosS(1);
    // pgtosS[0] = pgtoS;

    const int taskListSize = taskList.size();
#pragma omp parallel
    {
        int threadID = 0;
#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif // _OPENMP

#pragma omp for schedule(runtime)
        for (int i = 0; i < taskListSize; ++i) {
            const index_type shellIndexP = taskList[i].shellIndex1;
            const index_type shellIndexQ = taskList[i].shellIndex2;
            
            const int shellTypeP = orbitalInfo.getShellType(shellIndexP);
            const int shellTypeQ = orbitalInfo.getShellType(shellIndexQ);
            const int maxStepsP = 2 * shellTypeP + 1;
            const int maxStepsQ = 2 * shellTypeQ + 1;
            // const TlPosition posP = orbitalInfo.getPosition(shellIndexP);
            // const TlPosition posQ = orbitalInfo.getPosition(shellIndexQ);
            // const DfOverlapEngine::PGTOs pgtosP = DfOverlapEngine::getPGTOs(orbitalInfo, shellIndexP);
            // const DfOverlapEngine::PGTOs pgtosQ = DfOverlapEngine::getPGTOs(orbitalInfo, shellIndexQ);
            // const DfOverlapEngine::Query query(1, 0, 0, 0, shellTypeP, shellTypeQ, 0, 0);
            
            // this->pEngines_[threadID].calc0(query, posP, posQ, posR, posS, pgtosP, pgtosQ, pgtosR, pgtosS);
            this->pEngines_[threadID].calc(1, orbitalInfo, shellIndexP,
                                           0, orbitalInfo, shellIndexQ,
                                           0, orbitalInfo, -1,
                                           0, orbitalInfo, -1);
            
            int index = 0;
            // X
            for (int stepP = 0; stepP < maxStepsP; ++stepP) {
                const index_type globalShellIndexP = shellIndexP + stepP;

                for (int stepQ = 0; stepQ < maxStepsQ; ++stepQ) {
                    const index_type globalShellIndexQ = shellIndexQ + stepQ;

                    if ((shellIndexP != shellIndexQ) || (globalShellIndexP >= globalShellIndexQ)) {
                        const double dSdA = this->pEngines_[threadID].WORK[index];
                        const double dSdB = - dSdA;
                        
                        pMatX->add(globalShellIndexP, globalShellIndexQ, dSdA);
                        pMatX->add(globalShellIndexQ, globalShellIndexP, dSdB);
                    }
                    ++index;
                }
            }
            // Y
            for (int stepP = 0; stepP < maxStepsP; ++stepP) {
                const index_type globalShellIndexP = shellIndexP + stepP;

                for (int stepQ = 0; stepQ < maxStepsQ; ++stepQ) {
                    const index_type globalShellIndexQ = shellIndexQ + stepQ;

                    if ((shellIndexP != shellIndexQ) || (globalShellIndexP >= globalShellIndexQ)) {
                        const double dSdA = this->pEngines_[threadID].WORK[index];
                        const double dSdB = - dSdA;
                        
                        pMatY->add(globalShellIndexP, globalShellIndexQ, dSdA);
                        pMatY->add(globalShellIndexQ, globalShellIndexP, dSdB);
                    }
                    ++index;
                }
            }
            // Z
            for (int stepP = 0; stepP < maxStepsP; ++stepP) {
                const index_type globalShellIndexP = shellIndexP + stepP;

                for (int stepQ = 0; stepQ < maxStepsQ; ++stepQ) {
                    const index_type globalShellIndexQ = shellIndexQ + stepQ;

                    if ((shellIndexP != shellIndexQ) || (globalShellIndexP >= globalShellIndexQ)) {
                        const double dSdA = this->pEngines_[threadID].WORK[index];
                        const double dSdB = - dSdA;
                        
                        pMatZ->add(globalShellIndexP, globalShellIndexQ, dSdA);
                        pMatZ->add(globalShellIndexQ, globalShellIndexP, dSdB);
                    }
                    ++index;
                }
            }
        }
    }
}

