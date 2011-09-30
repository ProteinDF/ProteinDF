#include <queue>

#include "DfTwoElectronIntegral_Parallel.h"
#include "TlCommunicate.h"
#include "TlSparseSymmetricMatrix.h"
#include "TlFileSymmetricMatrix.h"
#include "TlFile.h"
#include "TlTime.h"

#define MASTER_TO_SLAVE_MSG_SIZE (3)

DfTwoElectronIntegral_Parallel::DfTwoElectronIntegral_Parallel(TlSerializeData* pPdfParam)
    : DfTwoElectronIntegral(pPdfParam)
{
    this->densityMatrixCacheMemSize_ = 100 * 1024UL * 1024UL; // 100MB
    this->isDebugOutMsg_ = (*pPdfParam)["debugout_TEI_parallel_msg"].getBoolean();
    
    this->isTeiOriginalCommunication_ = false;
    if ((*pPdfParam)["tei_original_communication"].getStr().empty() != true) {
        this->isTeiOriginalCommunication_ = (*pPdfParam)["tei_original_communication"].getBoolean();
    }
    
    this->numOfTasksPerProc_ = 2;
    if ((*pPdfParam)["xc_ms_job_per_proc"].getStr().empty() != true) {
        this->numOfTasksPerProc_ = (*pPdfParam)["xc_ms_job_per_proc"].getInt();
    }
}


DfTwoElectronIntegral_Parallel::~DfTwoElectronIntegral_Parallel()
{
//     std::cerr << "DfTwoElectronIntegral_Parallel::~DfTwoElectronIntegral_Parallel()" << std::endl;
}


void DfTwoElectronIntegral_Parallel::logger(const std::string& str) const
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    if (rComm.isMaster() == true) {
        DfTwoElectronIntegral::logger(str);
    }
}


void DfTwoElectronIntegral_Parallel::loggerP(const std::string& str) const
{
    DfTwoElectronIntegral::logger(str);
}


std::vector<DfTwoElectronIntegral::IKShellPair>
DfTwoElectronIntegral_Parallel::getLocalIKShellPairList(const int nIKShellPairIndex)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int nProc = rComm.getNumOfProc();
    const int nRank = rComm.getRank();

    const int nGlobalStart = 0;
    const int nGlobalEnd = this->m_IKShellPairList[nIKShellPairIndex].size();
    const int nRange = nGlobalEnd - nGlobalStart;
    const int nInterval = (nRange + (nProc -1)) / nProc; // +(nProc-1) は余り用
    const int nLocalStart = nGlobalStart + nRank * nInterval; // nProc = 0, 1, 2, ...
    const int nLocalEnd   = std::min((nGlobalStart + (nRank +1) * nInterval), nGlobalEnd);

    const std::vector<IKShellPair>::iterator pBegin = this->m_IKShellPairList[nIKShellPairIndex].begin();

    // for debug
//   {
//     std::cout << TlUtils::format("[%d] global(%d, %d) local(%d, %d)",
//               nRank, nGlobalStart, nGlobalEnd, nLocalStart, nLocalEnd) << std::endl;
//   }

//   std::vector<IKShellPair> aAnswer(pBegin +nLocalStart, pBegin +nLocalEnd);
    const int size = nLocalEnd - nLocalStart;
    std::vector<IKShellPair> aAnswer;
    if (size > 0) {
        aAnswer.resize(size);
        for (int i = 0; i < size; ++i) {
            aAnswer[i] = this->m_IKShellPairList[nIKShellPairIndex][nLocalStart +i];
        }
    }

    // for parallel version
    //return (this->m_IKShellPairList[nIKShellPairIndex]);
    return aAnswer;
}

void DfTwoElectronIntegral_Parallel::finalize(TlSymmetricMatrix& K)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    rComm.allReduce_SUM(K);
}


void DfTwoElectronIntegral_Parallel::getContractKMatrixByIntegralDriven(const TlSymmetricMatrix& P,
                                                                        TlSymmetricMatrix* pK)
{
    // initialize
    //const int nMaxShellType = 3; // s, p, d
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

    // Matrix
    const int nMaxIShell = this->m_TlOrbitalInfo.getNumOfOrbitals();
    pK->resize(nMaxIShell);

    //
    if (this->isMasterSlave_ == true) {
        this->loggerTime(" 2ERI calc. using Integral-Driven method(M/S).");
        this->getContractKMatrixByIntegralDrivenCore_MS(P, pK);
    } else {
        this->loggerTime(" 2ERI calc. using Integral-Driven method(D&C).");
        this->getContractKMatrixByIntegralDrivenCore(P, pK);
    }
}


void DfTwoElectronIntegral_Parallel::getContractKMatrixByIntegralDrivenCore_MS(const TlSymmetricMatrix& P,
                                                                               TlSymmetricMatrix* pK)
{
    assert(pK != NULL);

    TlCommunicate& rComm = TlCommunicate::getInstance();
    enum {
        REQUEST_JOB,
        ASSIGN_JOB,
        FINISH_JOB
    };

    if (rComm.isMaster() == true) {
        const int numOfProc = rComm.getNumOfProc();
        int numOfFinishJobProc = 0;

        // for job assign
        unsigned int eriType = 0;
        int iShell = 0;
        int kShell = 0;
        this->assignJob(&eriType, &iShell, &kShell, NULL, true);

        while (numOfFinishJobProc < (numOfProc -1)) {
            int msg = 0;
            int src = 0;
            rComm.receiveDataFromAnySource(msg, &src, TAG_TEI_REQUEST_TASK);
            assert(msg == REQUEST_JOB);

            if (this->assignJob(&eriType, &iShell, &kShell, NULL) == true) {
                const int assignJob = ASSIGN_JOB;
                rComm.sendData(assignJob, src, TAG_TEI_ASSIGN_TASK);

                std::vector<int> jobPacket(3);
                jobPacket[0] = static_cast<int>(eriType);
                jobPacket[1] = iShell;
                jobPacket[2] = kShell;
                rComm.sendData(jobPacket, src);
                if (this->isDebugOutMsg_ == true) {
                    this->loggerP(TlUtils::format(" [DEBUG] ERI:%d (i,k)=(%d,%d)\n", eriType, iShell, kShell));
                }
            } else {
                const int finishJob = FINISH_JOB;
                if (this->isDebugOutMsg_ == true) {
                    this->loggerP(TlUtils::format(" [DEBUG] [%d] finish (%d/%d)\n", rComm.getRank(), numOfFinishJobProc, rComm.getNumOfProc()));
                }
                rComm.sendData(finishJob, src);
                ++numOfFinishJobProc;
            }
        }
        
    } else {
        const int root = 0;
        const int requestJob = REQUEST_JOB;
        while (true) {
            rComm.sendData(requestJob, root, TAG_TEI_REQUEST_TASK);

            int assignJob = 0;
            rComm.receiveData(assignJob, root);
            if (assignJob == ASSIGN_JOB) {
                std::vector<int> jobPacket(3);
                rComm.receiveData(jobPacket, root, TAG_TEI_ASSIGN_TASK);
                const unsigned int eriType = static_cast<unsigned int>(jobPacket[0]);
                const int iShell = jobPacket[1];
                const int kShell = jobPacket[2];
                this->IntegralDrivenCore(eriType, iShell, kShell, P, pK);
            } else {
                if (this->isDebugOutMsg_ == true) {
                    this->loggerP(TlUtils::format(" [DEBUG] [%d] break\n", rComm.getRank(), rComm.getNumOfProc()));
                }
                break;
            }
        }
    }

    if (this->isDebugOutMsg_ == true) {
        this->loggerP(TlUtils::format(" [DEBUG] [%d] finalize\n", rComm.getRank()));
    }
    this->finalize(*pK);
}


bool DfTwoElectronIntegral_Parallel::assignJob(unsigned int* pEriType,
                                               int* pIShell, int* pKShell,
                                               double* pProgress,
                                               bool isInitialize)
{
    const int nMaxShellType = 3; // s, p, d
    static int iShellType = nMaxShellType -1;
    static int jShellType = nMaxShellType -1;
    static int kShellType = nMaxShellType -1;
    static int lShellType = nMaxShellType -1;

    static int lastIKShellPairIndex = -1;
    static std::vector<IKShellPair> ikShellPairList;
    static std::size_t ik = 0;
    static std::size_t total_IK_ShellPairs = 0; // 進捗計算用
    static std::size_t current_IK_ShellPairs = 0; // 進捗計算用

    if (isInitialize == true) {
        iShellType = nMaxShellType -1;
        jShellType = nMaxShellType -1;
        kShellType = nMaxShellType -1;
        lShellType = nMaxShellType -1;
        
        lastIKShellPairIndex = -1;
        ikShellPairList.clear();
        ik = 0;

        // 総数を計算
        total_IK_ShellPairs = 0;
        for (int i = 0; i < nMaxShellType; ++i) {
            for (int k = 0; k < nMaxShellType; ++k) {
                const int ikShellPairIndex = 4*i +k;
                total_IK_ShellPairs += DfTwoElectronIntegral::getLocalIKShellPairList(ikShellPairIndex).size();
            }
        }
        total_IK_ShellPairs *= 9;
        current_IK_ShellPairs = 0;
            
        
        return true;
    }
    
    // I: from (int IShell = 0) to (IShell < nMaxIShell)
    for (; iShellType >= 0; --iShellType) {
        //const std::vector<index_type> IShellList = this->m_ShellList[iShellType];

        // K: from (int nKShell = 0) to (nKShell <= nIShell)
        for (; kShellType >= 0; --kShellType) {

            const int ikShellPairIndex = 4*iShellType +kShellType;
            if (lastIKShellPairIndex != ikShellPairIndex) {
                ikShellPairList = DfTwoElectronIntegral::getLocalIKShellPairList(ikShellPairIndex);
                lastIKShellPairIndex = ikShellPairIndex;
            }

            // J: from (int nJShell = 0) to (nJShell <= nIShell)
            for (; jShellType >= 0; --jShellType) {

                // L: from (int nLShell = 0) to (nLShell <= nIShell)
                for (; lShellType >= 0; --lShellType) {

                    *pEriType = DfTEI::getEriType(iShellType, jShellType, kShellType, lShellType);

                    // i-Shell, k-Shell
                    const std::size_t maxIKShellPairIndex = ikShellPairList.size();
                    while (ik < maxIKShellPairIndex) {
                        *pIShell = ikShellPairList[ik].nIShell;
                        *pKShell = ikShellPairList[ik].nKShell;
                        ++ik;
                        ++current_IK_ShellPairs;
                        if (pProgress != NULL) {
                            *pProgress = ((double)current_IK_ShellPairs / (double)total_IK_ShellPairs) * 100.0;
                        }
                        
                        return true;
                    }
                    ik = 0;
                }
                lShellType = nMaxShellType -1;
            }
            jShellType = nMaxShellType -1;
        }
        kShellType = nMaxShellType -1;
    }

    return false;
}


// for RT
void DfTwoElectronIntegral_Parallel::getContractKMatrixByRTmethod(const TlSymmetricMatrix& P,
                                                                  TlSymmetricMatrix* pK)
{

    assert(pK != NULL);

    // initialize
    //const int nMaxShellType = 3; // s, p, d
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


    //if (this->isMasterSlave_ == true) {
    //this->getContractKMatrixByRTmethodCore_MS(P, pK);
    //} else {
    this->loggerTime(" 2ERI calc. using RT method(D&C).");
    this->getContractKMatrixByRTmethodCore(P, pK);
    //}
}


// for ScaLAPACK ===============================================================
void DfTwoElectronIntegral_Parallel::getContractKMatrixByRTmethod(const TlDistributeSymmetricMatrix& P,
                                                                  TlDistributeSymmetricMatrix* pK)
{
    //std::cerr << "DfTwoElectronIntegral_Parallel::getContractKMatrixByRTmethod()"<< std::endl;
    assert(pK != NULL);
    TlCommunicate& rComm = TlCommunicate::getInstance();
    assert(rComm.checkNonBlockingCommunications() == true);
    const int rank = rComm.getRank();
    const int numOfProcs = rComm.getNumOfProc();

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

    // 密度行列をベクトルに分解、各ノードに割り当てる
    // TODO: 高速化
    this->loggerTime(" distribute density matrix");

    const int numOfAOs = this->m_nNumOfAOs;
    this->logger(TlUtils::format(" density matrix cache: %lu byte\n", this->densityMatrixCacheMemSize_));
    TlSparseVectorMatrix tmpP(numOfAOs, numOfAOs, this->densityMatrixCacheMemSize_);
    std::map<int, int> densityVectorHoldingTable;
    {
        int proc = 0;
        for (int AO_index = 0; AO_index < numOfAOs; ) {
            const int shellType = this->m_TlOrbitalInfo.getShellType(AO_index);
            const int shells = shellType * 2 + 1;
            for (int shell = 0; shell < shells; ++shell) {
                const int AO = AO_index + shell;
                const TlVector v = P.getRowVector(AO);
                if (proc == rank) {
                    tmpP.setRowVector(AO, v);
                }
            }

//             if (this->isDebugOutMsg_ == true) {
//                 this->logger(TlUtils::format(" [DEBUG] %d-th density vector is stored in [%d]\n",
//                                              AO_index, proc));
//             }
            densityVectorHoldingTable[AO_index] = proc;
            AO_index += shells;
            proc = (proc +1) % numOfProcs;
        }
    }
    
    //std::cerr << "DfTwoElectronIntegral_Parallel::getContractKMatrixByRTmethod() part2"<< std::endl;
    this->loggerTime(" 2ERI calc. using RT method(ScaLAPACK).");

    TlSparseSymmetricMatrix sparseK(this->m_nNumOfAOs);
    if (this->isTeiOriginalCommunication_ == true) {
        // 短冊状行列を保持し、通信するタイプ
        if (rComm.isMaster() == true) {
            this->logger(" original communication.\n");
            this->RT_master(&tmpP, densityVectorHoldingTable);
        } else {
            this->RT_slave(&tmpP, &sparseK, densityVectorHoldingTable);
        }
    } else {
        // ブロック行列を通信するタイプ
        if (rComm.isMaster() == true) {
            this->logger(" common communication.\n");
            this->RT_master2(P);
        } else {
            this->RT_slave2(P, &sparseK);
        }
    }

    // finalize
    //std::cerr << "DfTwoElectronIntegral_Parallel::getContractKMatrixByRTmethod() finalize"<< std::endl;
    this->loggerTime(" merge partial matrix.");
    pK->mergeSparseMatrix(sparseK);

    this->loggerTime(" finish.");
    assert(rComm.checkNonBlockingCommunications() == true);
}


void DfTwoElectronIntegral_Parallel::RT_master(TlSparseVectorMatrix* pP,
                                               const std::map<int, int>& densityVectorHoldingTable)
{
    assert(pP != NULL);
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProc();
    
    // for job assign
    unsigned int eriType = 0;
    int iShell = 0;
    int kShell = 0;
    this->assignJob(&eriType, &iShell, &kShell, NULL, true); // "true" means intialize.

    // 送信すべきベクトルのリスト
    SendVectorTasksType sendVectorTasks;

    int cmdFromSlave[3];
    bool isWaitingFromSlaveCmd = false;
    int cmdFromSlave2[3];
    bool isWaitingFromSlaveCmd2 = false; // DensityVector用
    int cmdToSlave[3];
    bool isWaitingToSlaveCmd = false;
    int terminateJobProc = 0;
    bool isKeepLoop = true;

    while (isKeepLoop == true) {
        int src = 0;
        // phase1 (他プロセスからのコマンドを受ける) ---------------------------
        if (isWaitingFromSlaveCmd != true) {
            rComm.iReceiveDataFromAnySourceX(cmdFromSlave, 3, TAG_COMMAND);
            isWaitingFromSlaveCmd = true;
            if (this->isDebugOutMsg_ == true) {
                this->loggerP(" [DEBUG] [0] waiting command message from slave.\n");
            }
        }
        if (isWaitingFromSlaveCmd2 != true) {
            rComm.iReceiveDataFromAnySourceX(cmdFromSlave2, 3, TAG_REQUEST_VECTOR);
            isWaitingFromSlaveCmd2 = true;
            if (this->isDebugOutMsg_ == true) {
                this->loggerP(" [DEBUG] [0] waiting request vector message from slave.\n");
            }
        }

        // slaveからの要求
        if ((isWaitingFromSlaveCmd == true) &&
            (rComm.test(cmdFromSlave, &src) == true)) {
            rComm.wait(cmdFromSlave);
            isWaitingFromSlaveCmd = false;

            switch (cmdFromSlave[0]) {
            case REQUEST_JOB:
                if (this->assignJob(&eriType, &iShell, &kShell, NULL) == true) {
                    cmdToSlave[0] = (int)eriType;
                    cmdToSlave[1] = iShell;
                    cmdToSlave[2] = kShell;
                    if (this->isDebugOutMsg_ == true) {
                        this->loggerP(TlUtils::format(" [DEBUG] [0] assign job: (%d *|%d *)[%d] to [%d]\n",
                                                      iShell, kShell, eriType, src));
                    }

                    if (isWaitingToSlaveCmd == true) {
                        rComm.wait(cmdToSlave);
                        isWaitingToSlaveCmd = false;
                    }
                    rComm.iSendDataX(cmdToSlave, 3, src, TAG_COMMAND);
                    isWaitingToSlaveCmd = true;
                } else {
                    // finish job
                    ++terminateJobProc;
                    // メッセージは後で送る
                }
                break;
                
            default:
                std::cerr << "Program Error: " << __FILE__ << __LINE__ << std::endl;
                std::cerr << "unknown msg from " << src << std::endl;
                break;
            }
        }

        if (terminateJobProc >= numOfProcs -1) {
            // 全slaveノードに終了メッセージを送る
            std::vector<std::vector<int> > terminateMsgToSlave(numOfProcs -1);
            for (int slave = 0; slave < numOfProcs -1; ++slave) {
                if (this->isDebugOutMsg_ == true) {
                    this->loggerP(TlUtils::format(" [DEBUG] [0] send finish msg to [%d]\n", slave +1));
                }
                terminateMsgToSlave[slave].resize(3);
                terminateMsgToSlave[slave][0] = COMPLETION_JOB;
                rComm.iSendDataX(&(terminateMsgToSlave[slave][0]), 3, slave +1, TAG_COMMAND);
            }
            for (int slave = 0; slave < numOfProcs -1; ++slave) {
                rComm.wait(&(terminateMsgToSlave[slave][0]));
            }

            isKeepLoop = false; //break loop
        }

        // slaveからの要求(request vector)
        if ((isWaitingFromSlaveCmd2 == true) &&
            (rComm.test(cmdFromSlave2, &src) == true)) {
            rComm.wait(cmdFromSlave2);
            isWaitingFromSlaveCmd2 = false;

            // slaveからは密度行列ベクトルの要求が渡される
            assert(cmdFromSlave2[0] == REQUEST_DENS_VECTOR);
            if (this->isDebugOutMsg_ == true) {
                this->loggerP(TlUtils::format(" [DEBUG] [%d] receive request vector of %d-th from [%d].\n",
                                              rComm.getRank(),
                                              cmdFromSlave2[1],
                                              src));
            }
            sendVectorTasks.push_back(SendVectorInfo(src, cmdFromSlave2[1]));
        }
        
        // phase2 (ベクトル送信)
        this->sendDensityVector(*pP, &sendVectorTasks);
    }

    // 受信待ち状態の変数をキャンセル
    if (isWaitingToSlaveCmd == true) {
        rComm.cancel(cmdToSlave);
        isWaitingToSlaveCmd = false;
    }
    if (isWaitingFromSlaveCmd2 == true) {
        rComm.cancel(cmdFromSlave2);
        isWaitingFromSlaveCmd2 = false;
    }

    if (this->isDebugOutMsg_ == true) {
        this->logger(" [DEBUG] [0] finish RT_master()\n");
    }
}


void DfTwoElectronIntegral_Parallel::RT_slave(TlSparseVectorMatrix* pP,
                                              TlSparseSymmetricMatrix* pK,
                                              const std::map<int, int>& densityVectorHoldingTable)
{
    assert(pP != NULL);
    assert(pK != NULL);
    assert(pK->getNumOfRows() == this->m_nNumOfAOs);

    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int root = 0;
    
    CalcTasksType calcTasks;             // 計算すべきタスクのリスト
    SendVectorTasksType sendVectorTasks; // 送信すべきベクトルのリスト

    // 準備完了メッセージを送信
    int cmdToRoot[3];
    cmdToRoot[0] = REQUEST_JOB;
    rComm.sendDataX(cmdToRoot, 3, root, TAG_COMMAND);

    int cmdFromRoot[3];
    bool isWaitingFromRootCmd = false;  // rootからのメッセージを受信中(true)
    int cmdFromSlave[3];
    bool isWaitingFromSlaveCmd = false; // slaveからのメッセージを受信中(true)
    std::map<int, int*> requestDensVecCmdTbl; // AO番号:リクエストコマンドポインタ
    std::map<int, double*> recvDensVecCmdTbl; // AO番号:電子密度ベクトルポインタ
    bool isFinished = false;
    bool isLoopBreak = false;

#ifdef _OPENMP
    const int ompNestedBackup = omp_get_nested();
    omp_set_nested(1);
#endif // _OPENMP
    
#pragma omp parallel num_threads(2)
    {
#ifdef _OPENMP
        const int ompThreadId = omp_get_thread_num();
#else
        const int ompThreadId = 0;
#endif // _OPENMP
        
        if (ompThreadId == 0) {
            int src = 0;
            while (isLoopBreak == false) {
                // phase1 (他プロセスからのコマンドを受ける) ---------------------------
                if (isWaitingFromSlaveCmd != true) {
                    rComm.iReceiveDataFromAnySourceX(cmdFromSlave, 3, TAG_REQUEST_VECTOR);
                    isWaitingFromSlaveCmd = true;
                    if (this->isDebugOutMsg_ == true) {
                        this->loggerP(TlUtils::format(" [DEBUG] [%d] waiting message from slave.\n",
                                                      rComm.getRank()));
                    }
                }

                // slaveからの要求
                if ((isWaitingFromSlaveCmd == true) &&
                    (rComm.test(cmdFromSlave, &src) == true)) {
                    rComm.wait(cmdFromSlave);
                    isWaitingFromSlaveCmd = false;
                    
                    // slaveからは密度行列ベクトルの要求が渡される
                    assert(cmdFromSlave[0] == REQUEST_DENS_VECTOR);
                    if (this->isDebugOutMsg_ == true) {
                        this->loggerP(TlUtils::format(" [DEBUG] [%d] receive request vector of %d-th from [%d].\n",
                                                      rComm.getRank(),
                                                      cmdFromSlave[1],
                                                      src));
                    }
                    sendVectorTasks.push_back(SendVectorInfo(src, cmdFromSlave[1]));
                }

                // phase3 (ベクトル送信) -----------------------------------------------
                this->sendDensityVector(*pP, &sendVectorTasks);

#pragma omp flush(isLoopBreak)
                {
                }
            }
            
        } else {
            assert(ompThreadId == 1);

            while (isLoopBreak == false) {
                // int src = 0;
                // phase1 (他プロセスからのコマンドを受ける) ---------------------------
                if (isWaitingFromRootCmd != true) {
                    rComm.iReceiveDataX(cmdFromRoot, MASTER_TO_SLAVE_MSG_SIZE, root, TAG_COMMAND);
                    isWaitingFromRootCmd = true;
                    if (this->isDebugOutMsg_ == true) {
                        this->loggerP(TlUtils::format(" [DEBUG] [%d] waiting message from root.\n",
                                                      rComm.getRank()));
                    }
                }
                
                // rootからの要求
                if ((isWaitingFromRootCmd == true) &&
                    (rComm.test(cmdFromRoot) == true)) {
                    rComm.wait(cmdFromRoot);
                    isWaitingFromRootCmd = false;
                    
                    if (cmdFromRoot[0] == COMPLETION_JOB) {
                        // ループ終了メッセージ
                        isFinished = true;
                        if (this->isDebugOutMsg_ == true) {
                            this->loggerP(TlUtils::format(" [DEBUG] [%d] receive terminate message from 0.\n",
                                                          rComm.getRank()));
                        }
                    } else {
                        // masterからは計算すべきi,kシェルのインデックスが渡される
                        unsigned int eriType = (unsigned int)cmdFromRoot[0];
                        int i_shell = cmdFromRoot[1];
                        int k_shell = cmdFromRoot[2];
                        calcTasks.push_back(CalcShellPair(0,
                                                          this->m_TlOrbitalInfo,
                                                          this->m_nNumOfAOs,
                                                          eriType, i_shell, k_shell));
                        if (this->isDebugOutMsg_ == true) {
                            this->loggerP(TlUtils::format(" [DEBUG] [%d] push back calc job (%d * | %d *), %d)\n",
                                                          rComm.getRank(),
                                                          i_shell, k_shell, eriType));
                        }
                    }
                }

                // phase2 (2電子積分計算) ----------------------------------------------
                CalcTasksType::iterator ctEnd = calcTasks.end();
                for (CalcTasksType::iterator ct = calcTasks.begin(); ct != ctEnd; ) {
                    const unsigned int eriType = ct->eriType;
                    const int i_shell = ct->i_shell;
                    const int k_shell = ct->k_shell;
                    
                    // 保有していない電子密度があれば、他ノードに請求する
                    if (pP->findRow(i_shell) == false) {
                        this->requestDensityVector(densityVectorHoldingTable, i_shell,
                                                   &requestDensVecCmdTbl,
                                                   &recvDensVecCmdTbl);
                    }
                    if (pP->findRow(k_shell) == false) {
                        this->requestDensityVector(densityVectorHoldingTable, k_shell,
                                                   &requestDensVecCmdTbl,
                                                   &recvDensVecCmdTbl);
                    }
                    
                    // リクエストコマンドの送信が完了しているかチェック
                    this->cleanupRequestDensVecCmdTbl(&requestDensVecCmdTbl);

                    // 電子密度の受信が完了しているかチェック
                    this->cleanupRecvDensVecCmdTbl(&recvDensVecCmdTbl, pP);

                    // 準備が出来ていれば計算する
                    if ((pP->findRow(i_shell) == true) &&
                        (pP->findRow(k_shell) == true)) {
                        if (this->isDebugOutMsg_ == true) {
                            this->loggerP(TlUtils::format(" [DEBUG] [%d] calc (%d * | %d *), %d)\n",
                                                          rComm.getRank(),
                                                          i_shell, k_shell, eriType));
                        }
                        this->RT_Core(eriType, i_shell, k_shell, *pP, *pP, pK);
                        
                        // 計算終了メッセージをmaster(0)に送る
                        cmdToRoot[0] = REQUEST_JOB;
                        rComm.sendDataX(cmdToRoot, MASTER_TO_SLAVE_MSG_SIZE, root, TAG_COMMAND);
                        
                        ct = calcTasks.erase(ct);
                        continue;
                    }
                    ++ct;
                }
                
                // phase4 (終了判断) ---------------------------------------------------
                //if ((isFinished == true) && (calcTasks.empty() == true) && (sendVectorTasks.empty() == true)) {
                if (isFinished == true) {
                    //cmdToRoot[0] = TERMINATE_JOB;
                    //rComm.sendDataX(cmdToRoot, 3, root, TAG_COMMAND);
                    if (this->isDebugOutMsg_ == true) {
                        this->loggerP(TlUtils::format(" [DEBUG] [%d] send finish job message to root.\n",
                                                      rComm.getRank()));
                    }
                    isLoopBreak = true;
                }

#pragma omp flush(isLoopBreak)
                {
                }
            } // end while
        }
    } // end parallel

#ifdef _OPENMP
    omp_set_nested(ompNestedBackup);
#endif // _OPENMP
    

    if (isWaitingFromSlaveCmd == true) {
        rComm.cancel(cmdFromSlave);
        isWaitingFromSlaveCmd = false;
    }

    if (this->isDebugOutMsg_ == true) {
        this->loggerP(TlUtils::format(" [DEBUG] [%d] finish RT_slave()\n",
                                      rComm.getRank()));
    }
}


void DfTwoElectronIntegral_Parallel::RT_master2(const TlDistributeSymmetricMatrix& P)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    assert(rComm.checkNonBlockingCommunications() == true);

    const int numOfProcs = rComm.getNumOfProc();
    const int numOfTasksPerProc = this->numOfTasksPerProc_;
    
    // for job assign
    unsigned int eriType = 0;
    index_type iShell = 0;
    index_type kShell = 0;
    this->assignJob(&eriType, &iShell, &kShell, NULL, true); // "true" means intialize.

    // 送信すべきベクトルのリスト
    SendVectorTasksType sendVectorTasks;

    int cmdFromSlave;
    bool isWaitingFromSlaveCmd = false;

    int cmdToSlave[MASTER_TO_SLAVE_MSG_SIZE];
    bool isWaitingToSlaveCmd = false;

    int terminateJobProc = 0;
    int prevProgressStep = -1;
    bool isKeepLoop = true;
    while (isKeepLoop == true) {
        int src = 0;

        // 密度行列の配布処理
        P.getSparseMatrixX(NULL);
        
        // slaveからの要求
        if (isWaitingFromSlaveCmd != true) {
            rComm.iReceiveDataFromAnySource(cmdFromSlave, TAG_COMMAND);
            isWaitingFromSlaveCmd = true;
            if (this->isDebugOutMsg_ == true) {
                this->loggerP(" [DEBUG] [0] waiting command message from slave.\n");
            }
        }

        if ((isWaitingFromSlaveCmd == true) &&
            (rComm.test(cmdFromSlave, &src) == true)) {
            rComm.wait(cmdFromSlave);
            isWaitingFromSlaveCmd = false;

            double progress = 0.0;
            if (this->assignJob(&eriType, &iShell, &kShell, &progress) == true) {
                cmdToSlave[0] = (int)eriType;
                cmdToSlave[1] = iShell;
                cmdToSlave[2] = kShell;

                // progress
                const int progressStep = int(progress / 10.0);
                if (prevProgressStep < progressStep) {
                    this->loggerTime(TlUtils::format("  %3d%% done.", progressStep * 10));
                    prevProgressStep = progressStep;
                }

                // debug
                if (this->isDebugOutMsg_ == true) {
                    this->loggerP(TlUtils::format(" [DEBUG] [0] assign job: (%d *|%d *)[%d] to [%d]\n",
                                                  iShell, kShell, eriType, src));
                }
                
                if (isWaitingToSlaveCmd == true) {
                    rComm.wait(cmdToSlave);
                    isWaitingToSlaveCmd = false;
                }
                rComm.iSendDataX(cmdToSlave, MASTER_TO_SLAVE_MSG_SIZE, src, TAG_COMMAND);
                isWaitingToSlaveCmd = true;
            } else {
                // finish job
                ++terminateJobProc;
                // ループから抜けた直後に終了メッセージを送る
                if (terminateJobProc >= ((numOfProcs -1) * this->numOfTasksPerProc_)) {
                    isKeepLoop = false;
                    break;
                }
            }
        }
    }


    // 全slaveに終了メッセージを送る
    if (isWaitingToSlaveCmd == true) {
        rComm.wait(cmdToSlave);
        isWaitingToSlaveCmd = false;
    }
    cmdToSlave[0] = COMPLETION_JOB;
    for (int rank = 1; rank <= (numOfProcs-1); ++rank) {
        for (int task = 0; task < numOfTasksPerProc; ++task) {
            rComm.sendDataX(cmdToSlave, MASTER_TO_SLAVE_MSG_SIZE, rank, TAG_COMMAND);
            if (this->isDebugOutMsg_ == true) {
                this->logger(TlUtils::format(" [DEBUG] [0] send terminate msg to %d.\n",
                                             rank));
            }
        }
    }

    // 受信待ち状態の変数をキャンセル
    if (isWaitingToSlaveCmd == true) {
        rComm.cancel(cmdToSlave);
        isWaitingToSlaveCmd = false;
    }

    P.getSparseMatrixX(NULL, true); // finalize
    

    if (this->isDebugOutMsg_ == true) {
        this->logger(" [DEBUG] [0] finish RT_master()\n");
    }
    assert(rComm.checkNonBlockingCommunications() == true);
}


void DfTwoElectronIntegral_Parallel::RT_slave2(const TlDistributeSymmetricMatrix& P,
                                               TlSparseSymmetricMatrix* pK)
{
    assert(pK != NULL);
    assert(pK->getNumOfRows() == this->m_nNumOfAOs);

    TlCommunicate& rComm = TlCommunicate::getInstance();
    assert(rComm.checkNonBlockingCommunications() == true);
    const int root = 0;
    const int numOfTasksPerProc = this->numOfTasksPerProc_;
    
    CalcTasksType calcTasks;             // 計算すべきタスクのリスト
    SendVectorTasksType sendVectorTasks; // 送信すべきベクトルのリスト

    std::vector<int> cmdToRoot(numOfTasksPerProc_, REQUEST_JOB);
    std::vector<bool> isWaitingMsgToRoot(numOfTasksPerProc_, false);
    for (int i = 0; i < numOfTasksPerProc_; ++i) {
        rComm.sendData(cmdToRoot[i], root, TAG_COMMAND);    // 準備完了メッセージを送信
    }

    std::vector<std::vector<int> > cmdFromRoot(this->numOfTasksPerProc_);
    for (int i = 0; i < numOfTasksPerProc_; ++i) {
        cmdFromRoot[i].resize(MASTER_TO_SLAVE_MSG_SIZE);
    }
    std::vector<bool> isWaitingFromRootCmd(numOfTasksPerProc_, false);  // rootからのメッセージを受信中(true)
    
    std::vector<bool> isFinishTasks(numOfTasksPerProc_, false);
    bool isFinished = false;
    bool isBreakLoop = false;

    // statics
//     TlTime time_all;
//     time_all.start();
//     TlTime time_calc;
//     time_calc.stop();
    int numOfChances = 0;
    int numOfMisHits = 0;

#ifdef _OPENMP
    const int ompNestedBackup = omp_get_nested();
    omp_set_nested(1);
#endif // _OPENMP
    
    while (isBreakLoop == false) {
        P.getSparseMatrixX(NULL);
        
        // rootからの要求
        for (int i = 0; i < numOfTasksPerProc; ++i) {
            if ((isFinishTasks[i] == false) &&
                (isWaitingFromRootCmd[i] != true)) {
                rComm.iReceiveDataX(&(cmdFromRoot[i][0]),
                                    MASTER_TO_SLAVE_MSG_SIZE,
                                    root, TAG_COMMAND);
                isWaitingFromRootCmd[i] = true;
                if (this->isDebugOutMsg_ == true) {
                    this->loggerP(TlUtils::format(" [DEBUG] [%d] waiting message from root.\n",
                                                  rComm.getRank()));
                }
            }
        }
        
        for (int i = 0; i < numOfTasksPerProc; ++i) {
            if ((isWaitingFromRootCmd[i] == true) &&
                (rComm.test(&(cmdFromRoot[i][0])) == true)) {
                rComm.wait(&(cmdFromRoot[i][0]));
                isWaitingFromRootCmd[i] = false;
                
                if (cmdFromRoot[i][0] == COMPLETION_JOB) {
                    // ループ終了メッセージ
                    isFinishTasks[i] = true;
                    
                    if (this->isDebugOutMsg_ == true) {
                        this->loggerP(TlUtils::format(" [DEBUG] [%d] receive terminate message from 0.\n",
                                                      rComm.getRank()));
                    }
                } else {
                    // masterからは計算すべきi,kシェルのインデックスが渡される
                    unsigned int eriType = (unsigned int)cmdFromRoot[i][0];
                    int i_shell = cmdFromRoot[i][1];
                    int k_shell = cmdFromRoot[i][2];
                    calcTasks.push_back(CalcShellPair(i,
                                                      this->m_TlOrbitalInfo,
                                                      this->m_nNumOfAOs,
                                                      eriType, i_shell, k_shell));
                    if (this->isDebugOutMsg_ == true) {
                        this->loggerP(TlUtils::format(" [DEBUG] [%d] push back calc job (%d * | %d *), %d)\n",
                                                      rComm.getRank(),
                                                      i_shell, k_shell, eriType));
                    }
                }
            }
        }
        
        {
            // 2電子積分計算
            CalcTasksType::iterator ctEnd = calcTasks.end();
            for (CalcTasksType::iterator ct = calcTasks.begin(); ct != ctEnd; ) {
                const unsigned int eriType = ct->eriType;
                const index_type i_shell = ct->i_shell;
                const index_type k_shell = ct->k_shell;
                
                if (ct->isReadyForPi == false) {
                    ct->isReadyForPi = P.getSparseMatrixX(&(ct->Pi));
                    ++numOfChances;
                    if (ct->isReadyForPi != true) {
                        ++numOfMisHits;
                    } else {
                        if (this->isDebugOutMsg_ == true) {
                            this->loggerP(TlUtils::format(" [DEBUG] [%d] P_i (%d, *) OK!\n",
                                                          rComm.getRank(), i_shell));
                        }
                    }
                }
                if (ct->isReadyForPk == false) {
                    ct->isReadyForPk = P.getSparseMatrixX(&(ct->Pk));
                    ++numOfChances;
                    if (ct->isReadyForPk != true) {
                        ++numOfMisHits;
                    } else {
                        if (this->isDebugOutMsg_ == true) {
                            this->loggerP(TlUtils::format(" [DEBUG] [%d] P_k (%d, *) OK!\n",
                                                          rComm.getRank(), k_shell));
                        }
                    }
                }
                
                // 準備が出来ていれば計算する
                if ((ct->isReadyForPi == true) && (ct->isReadyForPk == true)) {
                    if (this->isDebugOutMsg_ == true) {
                        this->loggerP(TlUtils::format(" [DEBUG] [%d] calc (%d * | %d *), %d)\n",
                                                      rComm.getRank(),
                                                      i_shell, k_shell, eriType));
                    }
                    //time_calc.start();
                    this->RT_Core(eriType, i_shell, k_shell,
                                  ct->Pi, ct->Pk, pK);
                    //time_calc.stop();
                    
                    // 計算終了メッセージをmaster(0)に送る
                    const int sessionId = ct->sessionId;
                    cmdToRoot[sessionId] = REQUEST_JOB;
                    if (isWaitingMsgToRoot[sessionId] == true) {
                        rComm.wait(cmdToRoot[sessionId]);
                        isWaitingMsgToRoot[sessionId] = false;
                    }
                    rComm.iSendData(cmdToRoot[sessionId], root, TAG_COMMAND);
                    isWaitingMsgToRoot[sessionId] = true;
                    if (this->isDebugOutMsg_ == true) {
                        this->loggerP(TlUtils::format(" [DEBUG] [%d] send request msg to root.n",
                                                      rComm.getRank()));
                    }
                    
                    ct = calcTasks.erase(ct);
                    continue;
                }
                ++ct;
            }
        }
        
        {
            // 全タスク終了判断
            isFinished = true;
            for (int i = 0; i < numOfTasksPerProc; ++i) {
                if (isFinishTasks[i] != true) {
                    isFinished = false;
                    break;
                }
            }
            
            // 終了判断
            if ((isFinished == true) &&
                (calcTasks.empty() == true) &&
                (sendVectorTasks.empty() == true)) {
                if (this->isDebugOutMsg_ == true) {
                    this->loggerP(TlUtils::format(" [DEBUG] [%d] send finish job message to root.\n",
                                                  rComm.getRank()));
                }
                //break;
                isBreakLoop = true;
            }
        }
        
    } // end while

#ifdef _OPENMP
    omp_set_nested(ompNestedBackup);
#endif // _OPENMP

//     time_all.stop();
//     this->loggerP(TlUtils::format(" [%5d] %f sec. / %f sec. (%5.2f%%) distribution: %d / %d (%5.2f%%)\n",
//                                   rComm.getRank(),
//                                   time_calc.getElapseTime(), time_all.getElapseTime(),
//                                   time_calc.getElapseTime() / time_all.getElapseTime() * 100.0,
//                                   numOfMisHits, numOfChances, (double(numOfMisHits) / double(numOfChances)) * 100.0));
    
    for (int i = 0; i < numOfTasksPerProc_; ++i) {
        if (isWaitingFromRootCmd[i] == true) {
            rComm.cancel(&(cmdFromRoot[i][0]));
            isWaitingFromRootCmd[i] = false;
        }
        if (isWaitingMsgToRoot[i] == true) {
            rComm.wait(cmdToRoot[i]);
            isWaitingMsgToRoot[i] = false;
        }
    }

    P.getSparseMatrixX(NULL, true); // finalize
    
    if (this->isDebugOutMsg_ == true) {
        this->loggerP(TlUtils::format(" [DEBUG] [%d] finish RT_slave()\n",
                                      rComm.getRank()));
    }
    assert(rComm.checkNonBlockingCommunications() == true);
}


void DfTwoElectronIntegral_Parallel::requestDensityVector(const std::map<int, int>& densityVectorHoldingTable,
                                                          const int shellBegin,
                                                          std::map<int, int*>* pRequestDensVecCmdTbl,
                                                          std::map<int, double*>* pRecvDensVecCmdTbl)
{
    assert(pRequestDensVecCmdTbl != NULL);
    if ((pRequestDensVecCmdTbl->find(shellBegin) == pRequestDensVecCmdTbl->end()) &&
        (pRecvDensVecCmdTbl->find(shellBegin) == pRecvDensVecCmdTbl->end())) {
        TlCommunicate& rComm = TlCommunicate::getInstance();

        std::map<int, int>::const_iterator it = densityVectorHoldingTable.find(shellBegin);
        assert(it != densityVectorHoldingTable.end());
        const int holder = it->second;

        int* pCmdToSlave = new int[3];
        pCmdToSlave[0] = REQUEST_DENS_VECTOR;
        pCmdToSlave[1] = shellBegin;
        pCmdToSlave[2] = 0; // offset
        if (this->isDebugOutMsg_ == true) {
            this->loggerP(TlUtils::format(" [DEBUG] [%d] send request density vector of %d-th to [%d].\n",
                                          rComm.getRank(), shellBegin, holder));
        }
        rComm.iSendDataX(pCmdToSlave, 3, holder, TAG_REQUEST_VECTOR);
        (*pRequestDensVecCmdTbl)[shellBegin] = pCmdToSlave;

        const int numOfAOs = this->m_nNumOfAOs;
        const int shells = this->m_TlOrbitalInfo.getShellType(shellBegin) * 2 + 1;
        double* pDensVec = new double[numOfAOs * shells];
        rComm.iReceiveDataX(pDensVec, (numOfAOs * shells), holder, TAG_VECTOR_CONTENTS);
        (*pRecvDensVecCmdTbl)[shellBegin] = pDensVec;
    }
}


void DfTwoElectronIntegral_Parallel::cleanupRequestDensVecCmdTbl(std::map<int, int*>* pRequestDensVecCmdTbl)
{
    assert(pRequestDensVecCmdTbl != NULL);
    TlCommunicate& rComm = TlCommunicate::getInstance();

    std::map<int, int*>::iterator itEnd = pRequestDensVecCmdTbl->end();
    for (std::map<int, int*>::iterator it = pRequestDensVecCmdTbl->begin(); it != itEnd; ) {
        if (rComm.test(*(it->second)) == true) {
            rComm.wait(*(it->second));
            if (this->isDebugOutMsg_ == true) {
                this->loggerP(TlUtils::format(" [DEBUG] [%d] finish sending request density vector of %d-th.\n",
                                              rComm.getRank(), it->first));
            }

            delete[] it->second;
            it->second = NULL;

            pRequestDensVecCmdTbl->erase(it++);
            continue;
        }
        ++it;
    }
}


void DfTwoElectronIntegral_Parallel::cleanupRecvDensVecCmdTbl(std::map<int, double*>* pRecvDensVecCmdTbl,
                                                              TlSparseVectorMatrix* pP)
{
    assert(pRecvDensVecCmdTbl != NULL);
    TlCommunicate& rComm = TlCommunicate::getInstance();

    const int numOfAOs = this->m_nNumOfAOs;
    std::map<int, double*>::iterator itEnd = pRecvDensVecCmdTbl->end();
    for (std::map<int, double*>::iterator it = pRecvDensVecCmdTbl->begin(); it != itEnd; ) {
        if (rComm.test(*(it->second)) == true) {
            const int shellBegin = it->first;
            const int shells = this->m_TlOrbitalInfo.getShellType(shellBegin) * 2 + 1;
            rComm.wait(*(it->second));
            if (this->isDebugOutMsg_ == true) {
                this->loggerP(TlUtils::format(" [DEBUG] [%d] recv request density vector of %d-th.\n",
                                              rComm.getRank(), shellBegin));
            }

            for (int shell = 0; shell < shells; ++shell) {
                TlVector tmp(it->second + numOfAOs * shell, numOfAOs);
                pP->setRowVector(shellBegin + shell, tmp);
            }

            delete[] it->second;
            it->second = NULL;

            pRecvDensVecCmdTbl->erase(it++);
            continue;
        }
        ++it;
    }
}
    

void DfTwoElectronIntegral_Parallel::sendDensityVector(const TlSparseVectorMatrix& P,
                                                       SendVectorTasksType* pSendVectorTasks)
{
    assert(pSendVectorTasks != NULL);
    TlCommunicate& rComm = TlCommunicate::getInstance();
    
    SendVectorTasksType::iterator svEnd = pSendVectorTasks->end();
    for (SendVectorTasksType::iterator sv = pSendVectorTasks->begin(); sv != svEnd; ) {
        const int to = sv->to;
        const int shellBegin = sv->shell;
        if (sv->pBuf == NULL) {
            const int numOfAOs = this->m_nNumOfAOs;
            const int shells = this->m_TlOrbitalInfo.getShellType(shellBegin) * 2 + 1;
            sv->pBuf = new double[numOfAOs * shells];
            for (int shell = 0; shell < shells; ++shell) {
                const int base = shell * numOfAOs;
                const TlVector v = P.getRowVector(shellBegin + shell);
                for (int i = 0; i < numOfAOs; ++i) {
                    *(sv->pBuf + base + i) = v.get(i);
                }
            }
            if (this->isDebugOutMsg_ == true) {
                this->loggerP(TlUtils::format(" [DEBUG] [%d] send density vector(shell=%d) to [%d]\n",
                                              rComm.getRank(),
                                              shellBegin,
                                              to));
            }
            rComm.iSendDataX(sv->pBuf, numOfAOs * shells, to, TAG_VECTOR_CONTENTS);
        } else {
            if (rComm.test(*(sv->pBuf)) == true) {
                rComm.wait(*(sv->pBuf));

                delete[] sv->pBuf;
                sv->pBuf = NULL;
                sv = pSendVectorTasks->erase(sv);
                continue;
            }
        }

        ++sv;
    }
}


