#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

#include "DfEriX_Parallel.h"
#include "DfTaskCtrl_Parallel.h"
#include "TlCommunicate.h"
#include "TlSparseSymmetricMatrix.h"

DfEriX_Parallel::DfEriX_Parallel(TlSerializeData* pPdfParam) 
    : DfEriX(pPdfParam) {

    this->calcMode_ = (*pPdfParam)["ERI_calcmode"].getInt();
}


DfEriX_Parallel::~DfEriX_Parallel() 
{
}


void DfEriX_Parallel::logger(const std::string& str) const
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfEriX::logger(str);
    }
}

void DfEriX_Parallel::loggerTime_local(const std::string& str) const
{
    std::string out = str;
    int size = out.size();
    if (size > 0) {
        if (out[size -1] == '\n') {
            out.erase(size -1, 1);
        }

        const std::string timeStr = "[" + TlTime::getNow() + "]";
        TlUtils::pad(out, (72 - timeStr.length()), ' ');
        out += (timeStr + "\n");
        DfObject::logger(out);
    }
}

DfTaskCtrl* DfEriX_Parallel::getDfTaskCtrlObject() const
{
    DfTaskCtrl* pDfTaskCtrl = new DfTaskCtrl_Parallel(this->pPdfParam_);
    return pDfTaskCtrl;
}


void DfEriX_Parallel::finalize(TlMatrix* pMtx)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    rComm.allReduce_SUM(*pMtx);
}


void DfEriX_Parallel::finalize(TlSymmetricMatrix* pMtx)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    rComm.allReduce_SUM(*pMtx);
}


void DfEriX_Parallel::finalize(TlVector* pVct)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    rComm.allReduce_SUM(*pVct);
}


void DfEriX_Parallel::getJ_D(const TlVector& rho, TlDistributeSymmetricMatrix* pJ)
{
    assert(pJ != NULL);

    const TlOrbitalInfo orbitalInfo((*(this->pPdfParam_))["coordinates"],
                                    (*(this->pPdfParam_))["basis_set"]);
    const TlOrbitalInfo_Density orbitalInfo_Density((*(this->pPdfParam_))["coordinates"],
                                                    (*(this->pPdfParam_))["basis_set_auxD"]);
    const ShellArrayTable shellArrayTable = this->makeShellArrayTable(orbitalInfo);
    const ShellArrayTable shellArrayTable_Density = this->makeShellArrayTable(orbitalInfo_Density);

    const TlSparseSymmetricMatrix schwarzTable = this->makeSchwarzTable(orbitalInfo);

    pJ->resize(this->m_nNumOfAOs);
    TlSparseSymmetricMatrix tmpJ(this->m_nNumOfAOs);
    //pJ->zeroClear();
    
    DfTaskCtrl* pDfTaskCtrl = this->getDfTaskCtrlObject();
    std::vector<DfTaskCtrl::Task2> taskList;
    bool hasTask = pDfTaskCtrl->getQueue(orbitalInfo,
                                         false,
                                         this->grainSize_, &taskList, true);
    while (hasTask == true) {
        this->getJ_part(orbitalInfo,
                        orbitalInfo_Density,
                        shellArrayTable_Density,
                        taskList,
                        schwarzTable,
                        rho, &tmpJ);
        
        hasTask = pDfTaskCtrl->getQueue(orbitalInfo,
                                        false,
                                        this->grainSize_, &taskList);
    }
    
    //this->finalize(pJ);
    pDfTaskCtrl->cutoffReport();

    pJ->mergeSparseMatrix(tmpJ);
    
    delete pDfTaskCtrl;
    pDfTaskCtrl = NULL;
}


void DfEriX_Parallel::getJab_D(TlDistributeSymmetricMatrix* pJab)
{
    assert(pJab != NULL);
    const index_type numOfAuxDens = this->m_nNumOfAux;
    pJab->resize(numOfAuxDens);
    TlSparseSymmetricMatrix tmpJab(numOfAuxDens);
    
    const TlOrbitalInfo_Density orbitalInfo_Density((*(this->pPdfParam_))["coordinates"],
                                                    (*(this->pPdfParam_))["basis_set_auxD"]);

    const ShellArrayTable shellArrayTable = this->makeShellArrayTable(orbitalInfo_Density);

    DfTaskCtrl* pDfTaskCtrl = this->getDfTaskCtrlObject();
    std::vector<DfTaskCtrl::Task2> taskList;
    bool hasTask = pDfTaskCtrl->getQueue(orbitalInfo_Density,
                                         false,
                                         this->grainSize_, &taskList, true);
    while (hasTask == true) {
        this->getJab_part(orbitalInfo_Density,
                          taskList,
                          &tmpJab);

        hasTask = pDfTaskCtrl->getQueue(orbitalInfo_Density,
                                        false,
                                        this->grainSize_, &taskList);
    }

    //this->finalize(pJab);
    pDfTaskCtrl->cutoffReport();

    pJab->mergeSparseMatrix(tmpJab);
    
    delete pDfTaskCtrl;
    pDfTaskCtrl = NULL;
}


void DfEriX_Parallel::getJ_D(const TlDistributeSymmetricMatrix& P,
                             TlDistributeVector* pRho)
{
    switch (this->calcMode_) {
    case CalcMode_UsingLocalMatrix:
        this->getJ_D_local(P, pRho);
        break;
        
    default:
        this->getJ_D_BG(P, pRho);
        break;
    }
}

void DfEriX_Parallel::getJ_D_local(const TlDistributeSymmetricMatrix& P,
                                   TlDistributeVector* pRho)
{
    this->logger(" using local matrix for density matrix.\n");
    
    assert(pRho != NULL);
    //this->clearCutoffStats();
    TlCommunicate& rComm = TlCommunicate::getInstance();
    
    // カットオフ値の設定
    const double maxDeltaP = P.getMaxAbsoluteElement();
    if (maxDeltaP < 1.0) {
        this->cutoffThreshold_ /= std::fabs(maxDeltaP);
        this->logger(TlUtils::format(" new cutoff threshold = % e\n",
                                     this->cutoffThreshold_));
    }

    const TlOrbitalInfo orbitalInfo((*(this->pPdfParam_))["coordinates"],
                                    (*(this->pPdfParam_))["basis_set"]);
    const TlOrbitalInfo_Density orbitalInfo_Density((*(this->pPdfParam_))["coordinates"],
                                                    (*(this->pPdfParam_))["basis_set_auxD"]);
    const ShellArrayTable shellArrayTable_Density = this->makeShellArrayTable(orbitalInfo_Density);

    const TlSparseSymmetricMatrix schwarzTable = this->makeSchwarzTable(orbitalInfo);

    TlMatrix localP;
    std::vector<index_type> rowIndexes;
    std::vector<index_type> colIndexes;
    this->expandLocalDensityMatrix(P, orbitalInfo, &localP, &rowIndexes, &colIndexes);
    const index_type numOfRowIndexes = rowIndexes.size();
    const index_type numOfColIndexes = colIndexes.size();
    
    std::vector<DfTaskCtrl::Task2> taskList;
    DfTaskCtrl::Task2 task;
    for (index_type p = 0; p < numOfRowIndexes; ) {
        const index_type shellIndexP = rowIndexes[p];
        task.shellIndex1 = shellIndexP;
        for (index_type q = 0; q < numOfColIndexes; ) {
            const index_type shellIndexQ = colIndexes[q];
            if (shellIndexP >= shellIndexQ) {
                task.shellIndex2 = shellIndexQ;
                taskList.push_back(task);
            }
            q += orbitalInfo.getShellType(shellIndexQ) * 2 + 1;
        }
        p += orbitalInfo.getShellType(shellIndexP) * 2 + 1;
    }
    
    TlVector tmpRho(this->m_nNumOfAux);
    this->loggerTime(" ERI start");
    this->getJ_part2(orbitalInfo,
                     orbitalInfo_Density,
                     shellArrayTable_Density,
                     taskList,
                     schwarzTable,
                     TlDistributeMatrix(P), &tmpRho);
    this->loggerTime_local(TlUtils::format("  waiting...#%6d", rComm.getRank()));
    this->loggerTime(" ERI end");

    // finalize
    rComm.allReduce_SUM(tmpRho);
    *pRho = TlDistributeVector(tmpRho);
    this->loggerTime(" finished");
}

void DfEriX_Parallel::getJ_D_BG(const TlDistributeSymmetricMatrix& P,
                                TlDistributeVector* pRho)
{
    this->logger(" background transportation for density matrix.\n");
    
    assert(pRho != NULL);
    //this->clearCutoffStats();
    TlCommunicate& rComm = TlCommunicate::getInstance();
    
    // カットオフ値の設定
    const double maxDeltaP = P.getMaxAbsoluteElement();
    if (maxDeltaP < 1.0) {
        this->cutoffThreshold_ /= std::fabs(maxDeltaP);
        this->logger(TlUtils::format(" new cutoff threshold = % e\n",
                                     this->cutoffThreshold_));
    }

    const TlOrbitalInfo orbitalInfo((*(this->pPdfParam_))["coordinates"],
                                    (*(this->pPdfParam_))["basis_set"]);
    const TlOrbitalInfo_Density orbitalInfo_Density((*(this->pPdfParam_))["coordinates"],
                                                    (*(this->pPdfParam_))["basis_set_auxD"]);
    const ShellArrayTable shellArrayTable_Density = this->makeShellArrayTable(orbitalInfo_Density);

    const TlSparseSymmetricMatrix schwarzTable = this->makeSchwarzTable(orbitalInfo);

    TlSparseSymmetricMatrix tmpP(this->m_nNumOfAOs);
    bool isSetTempP = false;
    TlVector tmpRho(this->m_nNumOfAux);

    DfTaskCtrl* pDfTaskCtrl = this->getDfTaskCtrlObject();
    std::vector<DfTaskCtrl::Task2> taskList;
    bool hasTask = pDfTaskCtrl->getQueue(orbitalInfo,
                                         false,
                                         this->grainSize_, &taskList, true);
    while (hasTask == true) {
        if (isSetTempP != true) {
            const int numOfTasks = taskList.size();
            for (int task = 0; task < numOfTasks; ++task) {
                const index_type shellIndexP = taskList[task].shellIndex1;
                const index_type shellIndexQ = taskList[task].shellIndex2;
                const int shellTypeP = orbitalInfo.getShellType(shellIndexP);
                const int shellTypeQ = orbitalInfo.getShellType(shellIndexQ);
                const int maxStepsP = 2 * shellTypeP + 1;
                const int maxStepsQ = 2 * shellTypeQ + 1;
                for (int i = 0; i < maxStepsP; ++i) {
                    const index_type indexP = shellIndexP + i;
                    for (int j = 0; j < maxStepsQ; ++j) {
                        const index_type indexQ = shellIndexQ + j;
                        tmpP.set(indexP, indexQ, 0.0);
                    }
                }
            }
            isSetTempP = true;
        }

        if (P.getSparseMatrixX(&tmpP, false) == true) {
            this->getJ_part(orbitalInfo,
                            orbitalInfo_Density,
                            shellArrayTable_Density,
                            taskList,
                            schwarzTable,
                            tmpP, &tmpRho);
            
            tmpP.zeroClear();
            isSetTempP = false;

            hasTask = pDfTaskCtrl->getQueue(orbitalInfo,
                                            false,
                                            this->grainSize_, &taskList);
            // std::cerr << TlUtils::format("[%d] DfEriX_Parallel::getJ_D() task=%s",
            //                              rComm.getRank(),
            //                              ((hasTask == true) ? "T" : "F"))
            //           << std::endl;
        }

        P.getSparseMatrixX(NULL, false);
    }
    
    // std::cerr << TlUtils::format("[%d] DfEriX_Parallel::getJ_D() loop end",
    //                              rComm.getRank())
    //           << std::endl;

    // int allProcFinished = 0;
    // if (rComm.isMaster() == true) {
    //     const int numOfProcs = rComm.getNumOfProcs();
    //     for (int proc = 1; proc < numOfProcs; ++proc) {
    //         rComm.sendData(allProcFinished, proc, TAG_ALL_PROC_FINISHED);
    //     }
    // } else {
    //     rComm.iReceiveData(allProcFinished, 0, TAG_ALL_PROC_FINISHED);
    //     while (true) {
    //         P.getSparseMatrixX(NULL, false);
            
    //         if (rComm.test(&allProcFinished) == true) {
    //             rComm.wait(&allProcFinished);
    //             break;
    //         }
    //     }
    // }
    this->waitAnotherProcs(P);
    P.getSparseMatrixX(NULL, true);

    // std::cerr << TlUtils::format("[%d] DfEriX_Parallel::getJ_D() end",
    //                              rComm.getRank())
    //           << std::endl;
    
    pDfTaskCtrl->cutoffReport();
    delete pDfTaskCtrl;
    pDfTaskCtrl = NULL;

    // finalize
    //this->finalize(pRho);
    rComm.allReduce_SUM(tmpRho);
    *pRho = TlDistributeVector(tmpRho);
}


void DfEriX_Parallel::getJpq_D(const TlDistributeSymmetricMatrix& P,
                               TlDistributeSymmetricMatrix* pJ)
{
    assert(pJ != NULL);
    //TlCommunicate& rComm = TlCommunicate::getInstance();
    pJ->resize(this->m_nNumOfAOs);

    const TlOrbitalInfo orbitalInfo((*(this->pPdfParam_))["coordinates"],
                                    (*(this->pPdfParam_))["basis_set"]);

    const TlSparseSymmetricMatrix schwarzTable = this->makeSchwarzTable(orbitalInfo);

    TlSparseSymmetricMatrix tmpP(this->m_nNumOfAOs);
    bool isSetTempP = false;
    TlSparseSymmetricMatrix tmpJ(this->m_nNumOfAOs);

    DfTaskCtrl* pDfTaskCtrl = this->getDfTaskCtrlObject();
    std::vector<DfTaskCtrl::Task4> taskList;
    bool hasTask = pDfTaskCtrl->getQueue4(orbitalInfo,
                                          schwarzTable,
                                          this->grainSize_, &taskList, true);
    while (hasTask == true) {
        if (isSetTempP != true) {
            const int numOfTasks = taskList.size();
            for (int task = 0; task < numOfTasks; ++task) {
                const index_type shellIndexP = taskList[task].shellIndex1;
                const index_type shellIndexQ = taskList[task].shellIndex2;
                const index_type shellIndexR = taskList[task].shellIndex3;
                const index_type shellIndexS = taskList[task].shellIndex4;
                const int shellTypeP = orbitalInfo.getShellType(shellIndexP);
                const int shellTypeQ = orbitalInfo.getShellType(shellIndexQ);
                const int shellTypeR = orbitalInfo.getShellType(shellIndexR);
                const int shellTypeS = orbitalInfo.getShellType(shellIndexS);
                const int maxStepsP = 2 * shellTypeP + 1;
                const int maxStepsQ = 2 * shellTypeQ + 1;
                const int maxStepsR = 2 * shellTypeR + 1;
                const int maxStepsS = 2 * shellTypeS + 1;
                for (int i = 0; i < maxStepsP; ++i) {
                    const index_type indexP = shellIndexP + i;
                    for (int j = 0; j < maxStepsQ; ++j) {
                        const index_type indexQ = shellIndexQ + j;
                        tmpP.set(indexP, indexQ, 0.0);
                    }
                }
                for (int i = 0; i < maxStepsR; ++i) {
                    const index_type indexR = shellIndexR + i;
                    for (int j = 0; j < maxStepsS; ++j) {
                        const index_type indexS = shellIndexS + j;
                        tmpP.set(indexR, indexS, 0.0);
                    }
                }
            }
            isSetTempP = true;
        }

        if (P.getSparseMatrixX(&tmpP, false) == true) {
            this->getJ_integralDriven_part(orbitalInfo,
                                           taskList,
                                           tmpP, &tmpJ);
            tmpP.zeroClear();
            isSetTempP = false;
            
            hasTask = pDfTaskCtrl->getQueue4(orbitalInfo,
                                             schwarzTable,
                                             this->grainSize_, &taskList);
        }

        P.getSparseMatrixX(NULL, false);
    }
                    
    // int allProcFinished = 0;
    // if (rComm.isMaster() == true) {
    //     const int numOfProcs = rComm.getNumOfProcs();
    //     for (int proc = 1; proc < numOfProcs; ++proc) {
    //         rComm.sendData(allProcFinished, proc, TAG_ALL_PROC_FINISHED);
    //     }
    // } else {
    //     rComm.iReceiveData(allProcFinished, 0, TAG_ALL_PROC_FINISHED);
    //     while (true) {
    //         P.getSparseMatrixX(NULL, false);
            
    //         if (rComm.test(&allProcFinished) == true) {
    //             rComm.wait(&allProcFinished);
    //             break;
    //         }
    //     }
    // }
    this->waitAnotherProcs(P);
    P.getSparseMatrixX(NULL, true);

    //this->finalize(pJ);
    pDfTaskCtrl->cutoffReport();

    delete pDfTaskCtrl;
    pDfTaskCtrl = NULL;

    pJ->mergeSparseMatrix(tmpJ);
}


void DfEriX_Parallel::getK_D(const TlDistributeSymmetricMatrix& P,
                             TlDistributeSymmetricMatrix* pK)
{
    switch (this->calcMode_) {
    case CalcMode_UsingLocalMatrix:
        this->getK_D_local(P, pK);
        break;
        
    default:
        this->getK_D_BG(P, pK);
        break;
    }
}


void DfEriX_Parallel::getK_D_BG(const TlDistributeSymmetricMatrix& P,
                                TlDistributeSymmetricMatrix* pK)
{
    this->logger(" background transportation for density matrix.\n");

    assert(pK != NULL);
    //TlCommunicate& rComm = TlCommunicate::getInstance();
    
    const index_type numOfAOs = this->m_nNumOfAOs;
    pK->resize(numOfAOs);
    
    const TlOrbitalInfo orbitalInfo((*(this->pPdfParam_))["coordinates"],
                                    (*(this->pPdfParam_))["basis_set"]);
    
    const TlSparseSymmetricMatrix schwarzTable = this->makeSchwarzTable(orbitalInfo);

    TlSparseSymmetricMatrix tmpP(this->m_nNumOfAOs);
    bool isSetTempP = false;
    TlSparseSymmetricMatrix tmpK(this->m_nNumOfAOs);

    DfTaskCtrl* pDfTaskCtrl = this->getDfTaskCtrlObject();
    std::vector<DfTaskCtrl::Task4> taskList;
    bool hasTask = pDfTaskCtrl->getQueue4(orbitalInfo,
                                          schwarzTable,
                                          this->grainSize_, &taskList, true);
    while (hasTask == true) {
        if (isSetTempP != true) {
            const int numOfTasks = taskList.size();
            for (int task = 0; task < numOfTasks; ++task) {
                const index_type shellIndexP = taskList[task].shellIndex1;
                const index_type shellIndexQ = taskList[task].shellIndex2;
                const index_type shellIndexR = taskList[task].shellIndex3;
                const index_type shellIndexS = taskList[task].shellIndex4;
                const int shellTypeP = orbitalInfo.getShellType(shellIndexP);
                const int shellTypeQ = orbitalInfo.getShellType(shellIndexQ);
                const int shellTypeR = orbitalInfo.getShellType(shellIndexR);
                const int shellTypeS = orbitalInfo.getShellType(shellIndexS);
                const int maxStepsP = 2 * shellTypeP + 1;
                const int maxStepsQ = 2 * shellTypeQ + 1;
                const int maxStepsR = 2 * shellTypeR + 1;
                const int maxStepsS = 2 * shellTypeS + 1;
                for (int i = 0; i < maxStepsQ; ++i) {
                    const index_type indexQ = shellIndexQ + i;
                    for (int j = 0; j < maxStepsS; ++j) {
                        const index_type indexS = shellIndexS + j;
                        tmpP.set(indexQ, indexS, 0.0);
                    }
                    for (int j = 0; j < maxStepsR; ++j) {
                        const index_type indexR = shellIndexR + j;
                        tmpP.set(indexQ, indexR, 0.0);
                    }
                }
                for (int i = 0; i < maxStepsP; ++i) {
                    const index_type indexP = shellIndexP + i;
                    for (int j = 0; j < maxStepsS; ++j) {
                        const index_type indexS = shellIndexS + j;
                        tmpP.set(indexP, indexS, 0.0);
                    }
                    for (int j = 0; j < maxStepsR; ++j) {
                        const index_type indexR = shellIndexR + j;
                        tmpP.set(indexP, indexR, 0.0);
                    }
                }
            }
            isSetTempP = true;
        }

        if (P.getSparseMatrixX(&tmpP, false) == true) {
            this->getK_integralDriven_part(orbitalInfo,
                                           taskList,
                                           tmpP, &tmpK);
            tmpP.zeroClear();
            isSetTempP = false;

            hasTask = pDfTaskCtrl->getQueue4(orbitalInfo,
                                             schwarzTable,
                                             this->grainSize_, &taskList);
            // if (hasTask == false) {
            //     std::cerr << TlUtils::format("[%d] DfEriX_Parallel::getK_D() taskEnd",
            //                                  rComm.getRank())
            //               << std::endl;
            // }
        }

        P.getSparseMatrixX(NULL, false);
    }
    // std::cerr << TlUtils::format("[%d] DfEriX_Parallel::getK_D() waiting",
    //                              rComm.getRank())
    //           << std::endl;
    
    // waiting another proc
    this->waitAnotherProcs(P);
    P.getSparseMatrixX(NULL, true);

    //this->finalize(pK);
    pDfTaskCtrl->cutoffReport();

    delete pDfTaskCtrl;
    pDfTaskCtrl = NULL;

    pK->mergeSparseMatrix(tmpK);
}


void DfEriX_Parallel::getK_D_local(const TlDistributeSymmetricMatrix& P,
                                   TlDistributeSymmetricMatrix* pK)
{
    this->logger(" using local matrix for density matrix.\n");

    assert(pK != NULL);
    TlCommunicate& rComm = TlCommunicate::getInstance();
    
    const index_type numOfAOs = this->m_nNumOfAOs;
    pK->resize(numOfAOs);
    
    const TlOrbitalInfo orbitalInfo((*(this->pPdfParam_))["coordinates"],
                                    (*(this->pPdfParam_))["basis_set"]);
    
    const TlSparseSymmetricMatrix schwarzTable = this->makeSchwarzTable(orbitalInfo);

    TlMatrix localP;
    std::vector<index_type> rowIndexes;
    std::vector<index_type> colIndexes;
    this->expandLocalDensityMatrix(P, orbitalInfo, &localP, &rowIndexes, &colIndexes);
    //const index_type numOfRowIndexes = rowIndexes.size();
    //const index_type numOfColIndexes = colIndexes.size();
    
    //DfTaskCtrl* pDfTaskCtrl = this->getDfTaskCtrlObject();
    std::vector<DfTaskCtrl::Task4> taskList;

    this->loggerTime(" ERI start");
    TlSparseSymmetricMatrix tmpK(this->m_nNumOfAOs);
    const TlDistributeMatrix tmpP(P);
    DfTaskCtrl dfTaskCtrl(this->pPdfParam_);

    bool hasTask = dfTaskCtrl.getQueue4_K0(orbitalInfo,
                                          schwarzTable,
                                          tmpP,
                                          rowIndexes, colIndexes,
                                          this->grainSize_, &taskList, true);
    while (hasTask == true) {
        DfEriX::getK_integralDriven_part(orbitalInfo,
                                         taskList,
                                         tmpP, &tmpK);
        hasTask = dfTaskCtrl.getQueue4_K0(orbitalInfo,
                                         schwarzTable,
                                         tmpP,
                                         rowIndexes, colIndexes,
                                         this->grainSize_, &taskList);
    }

    // bool hasTask = dfTaskCtrl.getQueue4_K(orbitalInfo,
    //                                       schwarzTable,
    //                                       tmpP,
    //                                       rowIndexes, colIndexes,
    //                                       this->grainSize_, &taskList, true);
    // dfTaskCtrl.getQueue4_K(orbitalInfo,
    //                        schwarzTable,
    //                        tmpP,
    //                        rowIndexes, colIndexes,
    //                        this->grainSize_, &taskList);
    //     hasTask = dfTaskCtrl.getQueue4_K(orbitalInfo,
    //                                      schwarzTable,
    //                                      tmpP,
    //                                      rowIndexes, colIndexes,
    //                                      this->grainSize_, &taskList);
    
    this->loggerTime_local(TlUtils::format("  waiting...#%6d", rComm.getRank()));
    this->loggerTime(" ERI end");

    pK->mergeSparseMatrix(tmpK);
    this->loggerTime(" finished");
}


void DfEriX_Parallel::waitAnotherProcs(const TlDistributeSymmetricMatrix& P)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    
    int allProcFinished = 0;
    if (rComm.isMaster() == true) {
        const int numOfProcs = rComm.getNumOfProcs();
        // recv
        int recvMsgCount = 0;
        bool isWaitingMsg = false;
        int proc = 0;
        while (true) {
            P.getSparseMatrixX(NULL, false);
            if (isWaitingMsg != true) {
                rComm.iReceiveDataFromAnySource(allProcFinished, TAG_ALL_PROC_FINISHED);
                isWaitingMsg = true;
            } else {
                if (rComm.test(&allProcFinished, &proc) == true) {
                    rComm.wait(&allProcFinished);
                    ++recvMsgCount;
                    isWaitingMsg = false;
                }
            }
            if (recvMsgCount >= (numOfProcs -1)) {
                if (isWaitingMsg == true) {
                    rComm.cancel(&allProcFinished);
                }
                break;
            }
        }
        // send
        for (int proc = 1; proc < numOfProcs; ++proc) {
            rComm.sendData(allProcFinished, proc, TAG_ALL_PROC_FINISHED);
        }
    } else {
        // send
        rComm.sendData(allProcFinished, 0, TAG_ALL_PROC_FINISHED);
        // recv
        rComm.iReceiveData(allProcFinished, 0, TAG_ALL_PROC_FINISHED);
        while (true) {
            P.getSparseMatrixX(NULL, false);
            
            if (rComm.test(&allProcFinished) == true) {
                rComm.wait(&allProcFinished);
                break;
            }
        }
    }
    // std::cerr << TlUtils::format("[%d] DfEriX_Parallel::getK_D() end",
    //                              rComm.getRank())
    //           << std::endl;
}


void DfEriX_Parallel::expandLocalDensityMatrix(const TlDistributeSymmetricMatrix& P,
                                               const TlOrbitalInfo& orbInfo,
                                               TlMatrix* pLocalP,
                                               std::vector<index_type>* pRowIndexes,
                                               std::vector<index_type>* pColIndexes)
{
    assert(pLocalP != NULL);

    const TlMatrix refP = TlDistributeMatrix(P).getLocalMatrix();
    const std::vector<index_type> refRowIndexes = P.getRowIndexTable();
    const std::vector<index_type> refColIndexes = P.getColIndexTable();
    const index_type numOfRefRowIndexes = refRowIndexes.size();
    const index_type numOfRefColIndexes = refColIndexes.size();

    *pRowIndexes = this->getExpandIndexes(refRowIndexes, orbInfo);
    *pColIndexes = this->getExpandIndexes(refColIndexes, orbInfo);
    const index_type numOfLocalRows = pRowIndexes->size();
    const index_type numOfLocalCols = pColIndexes->size();

    pLocalP->resize(numOfLocalRows, numOfLocalCols);
    index_type ref_r = 0;
    index_type ref_c = 0;
    for (index_type r = 0; r < numOfLocalRows; ++r) {
        const index_type globalRow = (*pRowIndexes)[r];

        bool isFoundRow = false;
        while ((ref_r < numOfRefRowIndexes) && (refRowIndexes[ref_r] <= globalRow)) {
            if (refRowIndexes[ref_r] == globalRow) {
                isFoundRow = true;
                break;
            }
            ++ref_r;
        }

        if (isFoundRow == true) {
            for (index_type c = 0; c < numOfLocalCols; ++c) {
                index_type globalCol = (*pColIndexes)[c];

                bool isFoundCol = false;
                while ((ref_c < numOfRefColIndexes) && (refColIndexes[ref_c] <= globalCol)) {
                    if (refColIndexes[ref_c] == globalCol) {
                        isFoundCol = true;
                        break;
                    }
                    ++ref_c;
                }

                if (isFoundCol == true) {
                    const double value = refP.get(ref_r, ref_c);
                    pLocalP->set(r, c, value);
                }
            }
        }
        
    }
}

std::vector<DfObject::index_type>
DfEriX_Parallel::getExpandIndexes(const std::vector<index_type>& refIndexes,
                                  const TlOrbitalInfo& orbInfo)
{
    static const index_type BlockStartIndexTable[] = {
        0,                // s
        0, -1, -2,        // p
        0, -1, -2, -3, -4 // d
    };

    const index_type numOfRefIndexes = refIndexes.size();
    std::vector<index_type> newIndexes;
    newIndexes.reserve(numOfRefIndexes * 5); // 全てがd軌道(dzz)だと仮定して5倍した

    index_type prevBlockStartIndex = -1;
    for (index_type i = 0; i < numOfRefIndexes; ++i) {
        const index_type orb = refIndexes[i];
        const int basisType = orbInfo.getBasisType(orb);
        assert(basisType < (sizeof(BlockStartIndexTable)/sizeof(index_type)));
        const index_type blockStartIndex = orb + BlockStartIndexTable[basisType];
        if (prevBlockStartIndex != blockStartIndex) {
            const int blockSize = orbInfo.getShellType(orb) * 2 + 1;
            for (int j = 0; j < blockSize; ++j) {
                newIndexes.push_back(blockStartIndex + j);
            }
            prevBlockStartIndex = blockStartIndex;
        }
    }

    return newIndexes;
}


void DfEriX_Parallel::getJ_part2(const TlOrbitalInfo& orbitalInfo,
                                 const TlOrbitalInfo_Density& orbitalInfo_Density,
                                 const ShellArrayTable& shellArrayTable_Density,
                                 const std::vector<DfTaskCtrl::Task2>& taskList,
                                 const TlSparseSymmetricMatrix& schwarzTable,
                                 const TlDistributeMatrix& P, TlVector* pRho)
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
                
                    const DfEriEngine::CGTO_Pair RS = this->pEriEngines_[threadID].getCGTO_pair(orbitalInfo_Density,
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
                                const double P_pq = coef * P.getLocal(indexP, indexQ);
                                
                                for (int k = 0; k < maxStepsR; ++k) {
                                    const index_type indexR = shellIndexR + k;
                                    
                                    const double value = this->pEriEngines_[threadID].WORK[index];
                                    pRho->add(indexR, P_pq * value);
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
    }
}

// void DfEriX_Parallel::getK_integralDriven_part2(const TlOrbitalInfoObject& orbitalInfo,
//                                                 const std::vector<DfTaskCtrl::Task4>& taskList,
//                                                 const TlDistributeMatrix& P, TlMatrixObject* pK)
// {
//     const int taskListSize = taskList.size();
//     const double pairwisePGTO_cutoffThreshold = this->cutoffEpsilon3_;

// #pragma omp parallel
//     {
//         int threadID = 0;
// #ifdef _OPENMP
//         threadID = omp_get_thread_num();
// #endif // _OPENMP

//         this->pEriEngines_[threadID].setPrimitiveLevelThreshold(this->cutoffEpsilon3_);

// #pragma omp for schedule(runtime)
//         for (int i = 0; i < taskListSize; ++i) {
//             const index_type shellIndexP = taskList[i].shellIndex1;
//             const index_type shellIndexQ = taskList[i].shellIndex2;
//             const index_type shellIndexR = taskList[i].shellIndex3;
//             const index_type shellIndexS = taskList[i].shellIndex4;
//             const int shellTypeP = orbitalInfo.getShellType(shellIndexP);
//             const int shellTypeQ = orbitalInfo.getShellType(shellIndexQ);
//             const int shellTypeR = orbitalInfo.getShellType(shellIndexR);
//             const int shellTypeS = orbitalInfo.getShellType(shellIndexS);
//             const int maxStepsP = 2 * shellTypeP + 1;
//             const int maxStepsQ = 2 * shellTypeQ + 1;
//             const int maxStepsR = 2 * shellTypeR + 1;
//             const int maxStepsS = 2 * shellTypeS + 1;
                        
//             const DfEriEngine::CGTO_Pair PQ = this->pEriEngines_[threadID].getCGTO_pair(orbitalInfo,
//                                                                                         shellIndexP,
//                                                                                         shellIndexQ,
//                                                                                         pairwisePGTO_cutoffThreshold);
//             const DfEriEngine::CGTO_Pair RS = this->pEriEngines_[threadID].getCGTO_pair(orbitalInfo,
//                                                                                         shellIndexR,
//                                                                                         shellIndexS,
//                                                                                         pairwisePGTO_cutoffThreshold);
//             const DfEriEngine::Query queryPQ(0, 0, shellTypeP, shellTypeQ);
//             const DfEriEngine::Query queryRS(0, 0, shellTypeR, shellTypeS);

//             this->pEriEngines_[threadID].calc(queryPQ, queryRS, PQ, RS);
                        
//             this->storeK_integralDriven2(shellIndexP, maxStepsP,
//                                          shellIndexQ, maxStepsQ,
//                                          shellIndexR, maxStepsR,
//                                          shellIndexS, maxStepsS,
//                                          this->pEriEngines_[threadID], P, pK);
//         }
//     }
// }

// void DfEriX_Parallel::storeK_integralDriven2(const index_type shellIndexP, const int maxStepsP,
//                                              const index_type shellIndexQ, const int maxStepsQ,
//                                              const index_type shellIndexR, const int maxStepsR,
//                                              const index_type shellIndexS, const int maxStepsS,
//                                              const DfEriEngine& engine,
//                                              const TlDistributeMatrix& P,
//                                              TlMatrixObject* pK)
// {
//     assert(pK != NULL);

//     int index = 0;
//     for (int i = 0; i < maxStepsP; ++i) {
//         const index_type indexP = shellIndexP + i;
                                    
//         for (int j = 0; j < maxStepsQ; ++j) {
//             const index_type indexQ = shellIndexQ + j;
//             if (indexP < indexQ) {
//                 // bypass
//                 index += (maxStepsR * maxStepsS);
//                 continue;
//             }
                                        
//             for (int k = 0; k < maxStepsR; ++k) {
//                 const index_type indexR = shellIndexR + k;
//                 if (shellIndexQ == shellIndexS) {
//                     if (indexP < indexR) {
//                         index += maxStepsS;
//                         continue;
//                     }
//                 }
                                            
//                 for (int l = 0; l < maxStepsS; ++l) {
//                     const index_type indexS = shellIndexS + l;
//                     if (indexS > ((indexP == indexR) ? indexQ : indexR)) {
//                         // bypass
//                         ++index;
//                         continue;
//                     }
                                                
//                     const double value = engine.WORK[index];
                                                
//                     // Eq.1 : (indexP, indexR) <= (indexQ, indexS)
//                     const double coefEq1 = ((indexP == indexR) && (indexQ != indexS)) ? 2.0 : 1.0;
//                     pK->add(indexP, indexR, - coefEq1 * P.getLocal(indexQ, indexS) * value);
// #ifdef DEBUG_K
//                     this->IA_K_ID1_.countUp(indexP, indexR, indexQ, indexS, coefEq1);
// #endif // DEBUG_K
                    
//                     if (indexR != indexS) { // Eq.1 != Eq.2
//                         // Eq.2 : (indexP, indexS) <= (indexQ, indexR)
//                         const double coefEq2 = ((indexP == indexS) && (indexQ != indexR)) ? 2.0 : 1.0;
//                         pK->add(indexP, indexS, - coefEq2 * P.getLocal(indexQ, indexR) * value);
// #ifdef DEBUG_K
//                         this->IA_K_ID2_.countUp(indexP, indexS, indexQ, indexR, coefEq2);
// #endif // DEBUG_K
//                     }
                                                
//                     if (indexP != indexQ) { // (Eq.1, Eq.2) != (Eq.3, Eq.4)
//                         // Eq.3 : (indexQ, indexR) <= (indexP, indexS)
//                         if ((indexP != indexR) || (indexQ != indexS)) { // Eq.2 != Eq.3 : !((indexP, indexS) == (indexQ, indexR))
//                             const double coefEq3 = ((indexQ == indexR) && (indexP != indexS)) ? 2.0 : 1.0;
//                             pK->add(indexQ, indexR, - coefEq3 * P.getLocal(indexP, indexS) * value);
// #ifdef DEBUG_K
//                             this->IA_K_ID3_.countUp(indexQ, indexR, indexP, indexS, coefEq3);
// #endif // DEBUG_K
//                         }
                                                    
//                         if (indexR != indexS) { // Eq.3 != Eq.4
//                             // Eq.4 : (indexQ, indexS) <= (indexP, indexR)
//                             const double coefEq4 = ((indexQ == indexS) && (indexP != indexR)) ? 2.0 : 1.0;
//                             pK->add(indexQ, indexS, - coefEq4 * P.getLocal(indexP, indexR) * value);
// #ifdef DEBUG_K
//                             this->IA_K_ID4_.countUp(indexQ, indexS, indexP, indexR, coefEq4);
//                             // std::cerr << TlUtils::format("K_EQ4 <%2d %2d %2d %2d>: (%2d %2d %2d %2d), coef=%d, value=%d",
//                             //                              shellIndexP, shellIndexQ, shellIndexR, shellIndexS,
//                             //                              indexQ, indexS, indexP, indexR,
//                             //                              (int)coefEq4,
//                             //                              this->IA_K_ID4_.getCount(indexQ, indexS, indexP, indexR))
//                             //           << std::endl;
// #endif // DEBUG_K
//                         }
//                     }
//                     ++index;
//                 }
//             }
                                        
//         }
//     }
// }


TlSparseSymmetricMatrix DfEriX_Parallel::makeSchwarzTable(const TlOrbitalInfoObject& orbitalInfo)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProcs();
    const int rank = rComm.getRank();
    
    const index_type maxShellIndex = orbitalInfo.getNumOfOrbitals();
    TlSparseSymmetricMatrix schwarz(maxShellIndex);

    DfEriEngine engine;
    engine.setPrimitiveLevelThreshold(0.0);

    int counter = 0;
    for (index_type shellIndexP = 0; shellIndexP < maxShellIndex; ) {
        const int shellTypeP = orbitalInfo.getShellType(shellIndexP);
        const int maxStepsP = 2 * shellTypeP + 1;

        for (index_type shellIndexQ = 0; shellIndexQ < maxShellIndex; ) {
            const int shellTypeQ = orbitalInfo.getShellType(shellIndexQ);
            const int maxStepsQ = 2 * shellTypeQ + 1;
                
            if (counter == rank) {
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
            }    

            shellIndexQ += maxStepsQ;
            
            ++counter;
            if (counter >= numOfProcs) {
                counter = 0;
            }
        }
        shellIndexP += maxStepsP;
    }

    rComm.gatherToMaster(schwarz);
    rComm.broadcast(schwarz);
    
    return schwarz;
}
