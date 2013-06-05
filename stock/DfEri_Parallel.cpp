#include <cstdlib>
#include <list>

#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

#include "DfEri_Parallel.h"
#include "TlCommunicate.h"
#include "TlSparseMatrix.h"
#include "TlSparseSymmetricMatrix.h"
#include "TlDistributeVector.h"
#include "TlFileSymmetricMatrix.h"
#include "TlFile.h"
#include "TlSparseVectorMatrix.h"
#include "TlOrbitalInfo.h"
#include "TlTime.h"

DfEri_Parallel::DfEri_Parallel(TlSerializeData* pPdfParam) : DfEri(pPdfParam)
{
    this->densityMatrixCacheMemSize_ = 100 * 1024UL * 1024UL; // 100MB
    
    this->numOfTasksPerProc_ = 4;
    if ((*pPdfParam)["tasks_per_proc"].getStr().empty() != true) {
        this->numOfTasksPerProc_ = (*pPdfParam)["tasks_per_proc"].getInt();
    }
    
    // this->isDebugOut_getDeltaT_ = true;
}


DfEri_Parallel::~DfEri_Parallel()
{
}


void DfEri_Parallel::cutoffReport()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int root = 0;
    
    std::size_t stat[12];
    stat[ 0] = this->CcountSS;
    stat[ 1] = this->CcountPS;
    stat[ 2] = this->CcountPP;
    stat[ 3] = this->CcountDS;
    stat[ 4] = this->CcountDP;
    stat[ 5] = this->CcountDD;
    stat[ 6] = this->TcountSS;
    stat[ 7] = this->TcountPS;
    stat[ 8] = this->TcountPP;
    stat[ 9] = this->TcountDS;
    stat[10] = this->TcountDP;
    stat[11] = this->TcountDD;
    
    rComm.reduce_SUM(stat, 12, root);

    this->CcountSS = stat[ 0];
    this->CcountPS = stat[ 1];
    this->CcountPP = stat[ 2];
    this->CcountDS = stat[ 3];
    this->CcountDP = stat[ 4];
    this->CcountDD = stat[ 5];
    this->TcountSS = stat[ 6];
    this->TcountPS = stat[ 7];
    this->TcountPP = stat[ 8];
    this->TcountDS = stat[ 9];
    this->TcountDP = stat[10];
    this->TcountDD = stat[11];
    
    DfEri::cutoffReport();
}


void DfEri_Parallel::getDeltaT(const TlSymmetricMatrix& rDeltaPpq, TlVector* pDeltaT)
{
    //std::cerr << "DfEri_Parallel::getDeltaT()" << std::endl;

    // Set new cutvalue
    {
        const double MAXdeltaPpq = rDeltaPpq.getMaxAbsoluteElement();
        if ((MAXdeltaPpq > this->m_dCutValue) && (MAXdeltaPpq < 1.0)) {
            this->m_dCutValue /= MAXdeltaPpq;
            this->logger(TlUtils::format(" new cutvalue is %.2e\n", this->m_dCutValue));
        } else {
            this->m_dCutValue = this->m_dStoredCutValue;
        }
    }

    DfEri::getDeltaT_core(&rDeltaPpq, pDeltaT);

    TlCommunicate& rComm = TlCommunicate::getInstance();
    rComm.allReduce_SUM(*pDeltaT);
}


std::vector<DfObject::IJShellPair>
DfEri_Parallel::getQueue(const TlOrbitalInfo* pOrbitalInfo,
                         const std::vector<std::vector<std::size_t> >& shellList,
                         const int maxScore, const bool initialize)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int nProc = rComm.getNumOfProc();
    const int nRank = rComm.getRank();

    static unsigned int count = 0;
    if (initialize == true) {
        count = 0;
        return DfEri::getQueue(pOrbitalInfo, shellList,
                               maxScore, initialize);
    }

    std::vector<IJShellPair> ij = DfEri::getQueue(pOrbitalInfo, shellList,
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


void DfEri_Parallel::finalizeIntegral(TlSymmetricMatrix& rMatrix)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    rComm.allReduce_SUM(rMatrix);
}

void DfEri_Parallel::finalizeIntegral(TlVector& rVector)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    rComm.allReduce_SUM(rVector);
}


// Sab =================================================================
void DfEri_Parallel::getSab(TlSymmetricMatrix* pSab)
{
    assert(pSab != NULL);
    pSab->resize(this->m_nNumOfAux);

    this->getSab_core(pSab);

    this->loggerTime(" finalize");
    TlCommunicate& rComm = TlCommunicate::getInstance();
    rComm.allReduce_SUM(*pSab);
}


// =============================================================================
void DfEri_Parallel::getdeltaHpqA(const TlVector& deltaRho,
                                  TlSymmetricMatrix& deltaH)
{
    DfEri::getdeltaHpqA(deltaRho, deltaH);
}


void DfEri_Parallel::getdeltaHpqA(const TlVector& deltaRho,
                                  TlDistributeSymmetricMatrix& deltaH)
{
    this->getDeltaHpqA_rev2(deltaRho, deltaH);
}


void DfEri_Parallel::getdeltaHpqA(const TlDistributeVector& deltaRho,
                                  TlDistributeSymmetricMatrix& deltaH)
{
    this->loggerTime(" gather distribute vector.");
    const TlVector globalDeltaRho(deltaRho.getVector());
    this->getdeltaHpqA(globalDeltaRho, deltaH);
}


////////////////////////////////////////////////////////////////////////////////
// ScaLAPACK
//

void DfEri_Parallel::getSab(TlDistributeSymmetricMatrix* pSab)
{
    assert(pSab != NULL);
    pSab->resize(this->m_nNumOfAux);
//     const std::size_t needMem = this->m_nNumOfAux * (this->m_nNumOfAux + 1) / 2 * sizeof(double);
//     if ((this->isWorkOnDisk_ == true) ||
//         (this->procMaxMemSize_ < needMem)) {
//         this->logger(" <alpha|beta> is build on disk.\n");
//         TlMatrix::useMemManager(true);
//     } else {
//         this->logger(" <alpha|beta> is build on memory.\n");
//         TlMatrix::useMemManager(false);
//     }

    TlSparseSymmetricMatrix Sab_part(this->m_nNumOfAux);
    this->getSab_core(&Sab_part);

    this->loggerTime(" finalize");
    pSab->mergeSparseMatrix(Sab_part);
}


void DfEri_Parallel::getDeltaT(const TlDistributeSymmetricMatrix& rDeltaPpq,
                               TlDistributeVector* pDeltaT)
{
    this->logger(" calculation of deltaT using distributed memory algorithm.\n");

    // Set new cutvalue
    {
        const double MAXdeltaPpq = rDeltaPpq.getMaxAbsoluteElement();
        if ((MAXdeltaPpq > this->m_dCutValue) && (MAXdeltaPpq < 1.0)) {
            this->m_dCutValue /= MAXdeltaPpq;
            this->logger(TlUtils::format(" new cutvalue in is %.2e\n", this->m_dCutValue));
        } else {
            this->m_dCutValue = this->m_dStoredCutValue;
        }
    }

    this->getDeltaT_core(rDeltaPpq, pDeltaT);
}


void DfEri_Parallel::getDeltaT_core(const TlDistributeSymmetricMatrix& P,
                                    TlDistributeVector* pDeltaT)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    assert(rComm.checkNonBlockingCommunications() == true);
    const int numOfProcs = rComm.getNumOfProc();
    const int numOfTasksPerProc = this->numOfTasksPerProc_;
    
    this->initializeCounter();
    this->loggerTime(" sigma_pq(P_pq * <pq|alpha>) start");

    // Master-Slave
    this->loggerTime(" calc method: M/S and ScaLAPACK");
    const static int JOB_PROTOCOL_SIZE = 2; // index=0: shell type, index=1: shell index
    TlVector tmpDeltaT(pDeltaT->getSize());
    if (rComm.isMaster() == true) {
        // MASTER ==============================================================
        {
            int shellType;
            index_type shell;
            DfEri::getQueueX(this->shellList_, &shellType, &shell, true); // initialize only
        }
            
        // slaveからのjobのリクエスト
        int cmdRequestJob = 0;
        bool isWaitingCmdRequestJob = false;

        // ジョブキュー
        typedef std::list<int> AssignJobQueueType;
        AssignJobQueueType assignJobQueue;

        std::vector<bool> isWaitingJobList(numOfProcs);
        std::vector<index_type*> jobList(numOfProcs);
        for (int i = 0; i < numOfProcs; ++i) {
            jobList[i] = new int[JOB_PROTOCOL_SIZE];
        }

        int terminateJobProc = 0;
        bool isTerminate = false;
        while (true) {
            P.getSparseMatrixX(NULL);

            int src = 0;
            // phase1
            if (isWaitingCmdRequestJob != true) {
                rComm.iReceiveDataFromAnySource(cmdRequestJob, TAG_REQUEST_JOB);
                isWaitingCmdRequestJob = true;
            }

            // phase2
            if ((isWaitingCmdRequestJob == true) &&
                (rComm.test(cmdRequestJob, &src) == true)) {
                rComm.wait(cmdRequestJob);
                isWaitingCmdRequestJob = false;

                assignJobQueue.push_back(src);
            }

            // phase3
            AssignJobQueueType::iterator itEnd = assignJobQueue.end();
            for (AssignJobQueueType::iterator it = assignJobQueue.begin(); it != itEnd; ) {
                const int slave = *it;

                int shellType = 0;
                index_type shellIndex = 0;
                const bool isContinued = DfEri::getQueueX(this->shellList_, &shellType, &shellIndex);
                if (isContinued == true) {
                    // assign job
                    if (isWaitingJobList[slave] == true) {
                        rComm.wait(jobList[slave]);
                        isWaitingJobList[slave] = false;
                    }
                    jobList[slave][0] = shellType;
                    jobList[slave][1] = shellIndex;
                    rComm.iSendDataX(jobList[slave], JOB_PROTOCOL_SIZE, slave, TAG_ASSIGN_JOB);
                    isWaitingJobList[slave] = true;
                } else {
                    // finish job
                    ++terminateJobProc;
                    if (terminateJobProc >= (rComm.getNumOfProc() -1) * numOfTasksPerProc) {
                        int cmdTerminateSlave = -1;
                        for (int rank = 0; rank < (rComm.getNumOfProc() -1); ++rank) {
                            for (int task = 0; task < numOfTasksPerProc; ++task) {
                                rComm.sendData(cmdTerminateSlave, rank +1, TAG_TERMINATE_SLAVE);
                            }
                        }
                        isTerminate = true;
                    }
                }

                it = assignJobQueue.erase(it);
            }

            // terminate
            if ((isTerminate == true) &&
                (assignJobQueue.empty() == true)) {
                break;
            }
        }

        // 後始末
        for (int i = 0; i < numOfProcs; ++i) {
            if (isWaitingJobList[i] == true) {
                rComm.wait(jobList[i]); // 送信はcancelではない
                isWaitingJobList[i] = false;
            }
            
            delete[] jobList[i];
            jobList[i] = NULL;
        }

        if (isWaitingCmdRequestJob == true) {
            rComm.cancel(cmdRequestJob);
            isWaitingCmdRequestJob = false;
        }
    } else {
        // SLAVE ===============================================================
        const int root = 0;
        // static const int TERMINATE_TASK = 1;
        
        // ジョブのリスト
        typedef std::list<DT_Job> JobListType;
        JobListType jobList;

        // タスク要求
        std::vector<int> cmdRequestTask(numOfTasksPerProc, 0);
        std::vector<bool> isWaitingCmdRequestTask(numOfTasksPerProc, false);
        for (int i = 0; i < numOfTasksPerProc; ++i) {
            rComm.iSendData(cmdRequestTask[i], root, TAG_REQUEST_JOB);
            isWaitingCmdRequestTask[i] = true;
        }
        
        // masterからのタスクの割り当て
        std::vector<std::vector<int> > cmdAssignTask(numOfTasksPerProc);
        std::vector<bool> isWaitingCmdAssignTask(numOfTasksPerProc);
        for (int i = 0; i < numOfTasksPerProc; ++i) {
            cmdAssignTask[i].resize(JOB_PROTOCOL_SIZE);
            isWaitingCmdAssignTask[i] = false;
        }
        
        // masterからの終了メッセージ
        std::vector<int> cmdTerminateSlave(numOfTasksPerProc, 0);
        std::vector<bool> isWaitingCmdTerminateSlave(numOfTasksPerProc, false);
        int terminateMsg = 0;

        bool isTerminate = false;

#ifdef _OPENMP
        // const int ompNestedBackup = omp_get_nested();
        omp_set_nested(1);
#endif // _OPENMP

        while (isTerminate == false) {
            P.getSparseMatrixX(NULL);
            
            // phase1
            for (int i = 0; i < numOfTasksPerProc; ++i) {
                if (isWaitingCmdAssignTask[i] != true) {
                    rComm.iReceiveDataX(&(cmdAssignTask[i][0]), JOB_PROTOCOL_SIZE, root, TAG_ASSIGN_JOB);
                    isWaitingCmdAssignTask[i] = true;
                }
            }
            for (int i = 0; i < numOfTasksPerProc; ++i) {
                if (isWaitingCmdTerminateSlave[i] != true) {
                    rComm.iReceiveData(cmdTerminateSlave[i], root, TAG_TERMINATE_SLAVE);
                    isWaitingCmdTerminateSlave[i] = true;
                }
            }
            
            // phase2
            for (int i = 0; i < numOfTasksPerProc; ++i) {
                if ((isWaitingCmdAssignTask[i] == true) &&
                    (rComm.test(&(cmdAssignTask[i][0])) == true)) {
                    rComm.wait(&(cmdAssignTask[i][0]));
                    isWaitingCmdAssignTask[i] = false;
                    
                    const int shellType = cmdAssignTask[i][0];
                    const int iShell = cmdAssignTask[i][1];
                    std::div_t d = std::div(shellType, DfEri::MAX_TYPE);
                    const int iShellType = d.quot;
                    const int jShellType = d.rem;
                    
                    TlSparseSymmetricMatrix partP(this->m_nNumOfAOs);
                    const int numOfJShells = this->shellList_[jShellType].size();
                    std::vector<IJShellPair> ijShellPairs;
                    ijShellPairs.reserve(numOfJShells);
                    for (int j = 0; j < numOfJShells; ++j) {
                        const index_type jShell = this->shellList_[jShellType][j];
                        
                        if ((iShellType == jShellType) && (iShell < jShell)) {
                            continue;
                        }
                        this->countupTotal(iShellType, jShellType);
                        if (this->isCutoffUsingSchwartzInequality(*(this->pOrbitalInfo_),
                                                                  iShell, jShell,
                                                                  this->m_dCutValue) == true) {
                            this->countupCutoff(iShellType, jShellType);
                            continue;
                        }
                        
                        ijShellPairs.push_back(IJShellPair(iShell, jShell));
                        
                        const int maxIOrbType = this->pOrbitalInfo_->getShellType(iShell) * 2 + 1;
                        const int maxJOrbType = this->pOrbitalInfo_->getShellType(jShell) * 2 + 1;
                        for (int iOrbType = 0; iOrbType < maxIOrbType; ++iOrbType) {
                            for (int jOrbType = 0; jOrbType < maxJOrbType; ++jOrbType) {
                                partP.set(iShell + iOrbType, jShell + jOrbType, 0.0);
                            }
                        }
                    }
                    
                    jobList.push_back(DT_Job(partP, ijShellPairs));
                }
            }
            for (int i = 0; i < numOfTasksPerProc; ++i) {
                if ((isWaitingCmdTerminateSlave[i] == true) &&
                    (rComm.test(cmdTerminateSlave[i]) == true)) {
                    rComm.wait(cmdTerminateSlave[i]);
                    isWaitingCmdTerminateSlave[i] = false;
                    
                    ++terminateMsg;
                }
            }
            
            // phase3
            JobListType::iterator itEnd = jobList.end();
            for (JobListType::iterator it = jobList.begin(); it != itEnd; ) {
                if (P.getSparseMatrixX(&(it->partP)) == true) {
                    // 計算
                    this->ericalcDT(it->ijShellPairs, &(it->partP), &tmpDeltaT);
                    
                    // 空いている通信ソケットをさがす
                    int socket = 0;
                    while (true) {
                        if (rComm.test(cmdRequestTask[socket]) == true) {
                            rComm.wait(cmdRequestTask[socket]);
                            isWaitingCmdRequestTask[socket] = false;
                            break;
                        }
                        ++socket;
                        if (socket > numOfTasksPerProc) {
                            socket = 0;
                        }
                    }
                    rComm.iSendData(cmdRequestTask[socket], root, TAG_REQUEST_JOB);
                    isWaitingCmdRequestTask[socket] = true;
                    
                    it = jobList.erase(it);
                }
                ++it;
            }
            
            if (jobList.empty() == true) {
                if (terminateMsg >= numOfTasksPerProc) {
                    isTerminate = true;
                    
                }
            }
        } // end while

        // 後始末
        for (int i = 0; i < numOfTasksPerProc; ++i) {
            if (isWaitingCmdAssignTask[i] == true) {
                rComm.cancel(&(cmdAssignTask[i][0]));
                isWaitingCmdAssignTask[i] = false;
            }
            
            if (isWaitingCmdTerminateSlave[i] == true) {
                rComm.cancel(cmdTerminateSlave[i]);
                isWaitingCmdTerminateSlave[i] = false;
            }

            if (isWaitingCmdRequestTask[i] == true) {
                rComm.cancel(cmdRequestTask[i]);
                isWaitingCmdRequestTask[i] = false;
            }
        }
    } // end else

    // 
    this->loggerTime(" finalize");
    P.getSparseMatrixX(NULL, true);
    rComm.allReduce_SUM(tmpDeltaT);
    *pDeltaT = tmpDeltaT;

    // cutoff report
    this->cutoffReport();

    this->loggerTime(" end");
    assert(rComm.checkNonBlockingCommunications() == true);
}


void DfEri_Parallel::getDeltaHpqA_rev2(const TlVector& deltaRho,
                                       TlDistributeSymmetricMatrix& deltaH)
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    assert(rComm.checkNonBlockingCommunications() == true);
    const int numOfProcs = rComm.getNumOfProc();

    this->initializeCounter();
    this->loggerTime(" J_pq = sigma_alpha(rho_alpha * <pq|alpha>): start");

    // Set new cutvalue
    const double MAXdeltaRho = deltaRho.getMaxAbsoluteElement();
    this->m_dCutValue = this->m_dStoredCutValue;
    if ((MAXdeltaRho > this->m_dCutValue) && (MAXdeltaRho < 1.0)) {
        this->m_dCutValue /= MAXdeltaRho;
        this->logger(TlUtils::format(" new cutvalue is %.2e\n", this->m_dCutValue));
    }

    this->loggerTime(" calc method: M/S and ScaLAPACK");
    const static int JOB_PROTOCOL_SIZE = 2; // index=0: shell type, index=1: shell index
    TlSparseSymmetricMatrix deltaHTmp(deltaH.getNumOfRows());
    if (rComm.isMaster() == true) {
        // MASTER ==============================================================
        {
            int shellType;
            index_type shell;
            DfEri::getQueueX(this->shellList_, &shellType, &shell, true); // initialize only
        }

        // slaveからのjobのリクエスト
        int cmdRequestJob = 0;
        bool isWaitingCmdRequestJob = false;

        // ジョブキュー
        typedef std::list<int> AssignJobQueueType;
        AssignJobQueueType assignJobQueue;

        //
        std::vector<bool> isWaitingJobList(numOfProcs);
        std::vector<index_type*> jobList(numOfProcs);
        for (int i = 0; i < numOfProcs; ++i) {
            jobList[i] = new int[JOB_PROTOCOL_SIZE];
        }

        // masterからslaveへの終了メッセージ
        int cmdTerminateSlave = 0;
        bool isWaitingCmdTerminateSlave = false;

        // slaveからの終了メッセージ
        int cmdTerminateOK = 0;
        bool isWaitingCmdTerminateOK = false;

        // 
        int terminateJobProc = 0;
        bool isTerminate = false;
        
        while (true) {
            int src = 0;
            // phase1
            if (isWaitingCmdRequestJob != true) {
                rComm.iReceiveDataFromAnySource(cmdRequestJob, TAG_REQUEST_JOB);
                isWaitingCmdRequestJob = true;
            }
            if (isWaitingCmdTerminateOK != true) {
                rComm.iReceiveDataFromAnySource(cmdTerminateOK, TAG_TERMINATE_OK);
                isWaitingCmdTerminateOK = true;
            }

            // phase2
            if ((isWaitingCmdRequestJob == true) &&
                (rComm.test(cmdRequestJob, &src) == true)) {
                rComm.wait(cmdRequestJob);
                isWaitingCmdRequestJob = false;

                assignJobQueue.push_back(src);
            }
            if ((isWaitingCmdTerminateOK == true) &&
                (rComm.test(cmdTerminateOK, &src) == true)) {
                rComm.wait(cmdTerminateOK);
                isWaitingCmdTerminateOK = false;
                ++terminateJobProc;
                if (terminateJobProc >= numOfProcs -1) {
                    isTerminate = true;
                }
            }

            // phase3
            AssignJobQueueType::iterator itEnd = assignJobQueue.end();
            for (AssignJobQueueType::iterator it = assignJobQueue.begin(); it != itEnd; ) {
                const int slave = *it;

                int shellType = 0;
                index_type shellIndex = 0;
                const bool isContinued = DfEri::getQueueX(this->shellList_, &shellType, &shellIndex);
                if (isContinued == true) {
                    // assign job
                    if (isWaitingJobList[slave] == true) {
                        rComm.wait(jobList[slave]);
                        isWaitingJobList[slave] = false;
                    }
                    jobList[slave][0] = shellType;
                    jobList[slave][1] = shellIndex;
                    rComm.iSendDataX(jobList[slave], JOB_PROTOCOL_SIZE, slave, TAG_ASSIGN_JOB);
                    isWaitingJobList[slave] = true;
                } else {
                    // finish job
                    if (isWaitingCmdTerminateSlave == true) {
                        rComm.wait(cmdTerminateSlave);
                        isWaitingCmdTerminateSlave = false;
                    }
                    rComm.iSendData(cmdTerminateSlave, slave, TAG_TERMINATE_SLAVE);
                    isWaitingCmdTerminateSlave = true;
                }

                it = assignJobQueue.erase(it);
            }

            // terminate
            if ((isTerminate == true) &&
                (assignJobQueue.empty() == true)) {
                break;
            }
        }

        // 後始末
        for (int i = 0; i < numOfProcs; ++i) {
            if (isWaitingJobList[i] == true) {
                rComm.wait(jobList[i]); // 送信はcancelではない
                isWaitingJobList[i] = false;
            }
            
            delete[] jobList[i];
            jobList[i] = NULL;
        }

        if (isWaitingCmdRequestJob == true) {
            rComm.cancel(cmdRequestJob);
            isWaitingCmdRequestJob = false;
        }
       
        if (isWaitingCmdTerminateSlave == true) {
            rComm.wait(cmdTerminateSlave); // 送信はcalcelではない
            isWaitingCmdTerminateSlave = false;
        }
        
        if (isWaitingCmdTerminateOK == true) {
            rComm.cancel(cmdTerminateOK);
            isWaitingCmdTerminateOK = false;
        }
    } else {
        // SLAVE ===============================================================
        const int root = 0;
        // masterからのジョブの割り当て
        int* pCmdAssignJob = new int[JOB_PROTOCOL_SIZE];
        bool isWaitingCmdAssignJob = false;
        // masterからの終了メッセージ
        int cmdTerminateSlave = 0;
        bool isWaitingCmdTerminateSlave = false;

        // ジョブのリスト
        typedef std::list<std::vector<IJShellPair> > JobListType;
        JobListType jobList;

        // 初期メッセージ
        int cmdRequestJob = 0;
        rComm.iSendData(cmdRequestJob, root, TAG_REQUEST_JOB);
        bool isWaitingCmdRequestJob = true;
        
        bool isTerminate = false;
        while (true) {
            // phase1
            if (isWaitingCmdAssignJob != true) {
                rComm.iReceiveDataX(pCmdAssignJob, JOB_PROTOCOL_SIZE, root, TAG_ASSIGN_JOB);
                isWaitingCmdAssignJob = true;
            }
            if (isWaitingCmdTerminateSlave != true) {
                rComm.iReceiveData(cmdTerminateSlave, root, TAG_TERMINATE_SLAVE);
                isWaitingCmdTerminateSlave = true;
            }

            // phase2
            if ((isWaitingCmdAssignJob == true) &&
                (rComm.test(pCmdAssignJob) == true)) {
                rComm.wait(pCmdAssignJob);
                isWaitingCmdAssignJob = false;

                const int shellType = pCmdAssignJob[0];
                const int iShell = pCmdAssignJob[1];
                std::div_t d = std::div(shellType, DfEri::MAX_TYPE);
                const int iShellType = d.quot;
                const int jShellType = d.rem;
                
                const int numOfJShells = this->shellList_[jShellType].size();
                std::vector<IJShellPair> ijShellPairs;
                ijShellPairs.reserve(numOfJShells);
                for (int j = 0; j < numOfJShells; ++j) {
                    const index_type jShell = this->shellList_[jShellType][j];

                    if ((iShellType == jShellType) && (iShell < jShell)) {
                        continue;
                    }
                    this->countupTotal(iShellType, jShellType);
                    if (this->isCutoffUsingSchwartzInequality(*(this->pOrbitalInfo_),
                                                              iShell, jShell,
                                                              this->m_dCutValue) == true) {
                        this->countupCutoff(iShellType, jShellType);
                        continue;
                    }
                    ijShellPairs.push_back(IJShellPair(iShell, jShell));
                }
                
                jobList.push_back(ijShellPairs);
            }
            if ((isWaitingCmdTerminateSlave == true) &&
                (rComm.test(cmdTerminateSlave) == true)) {
                rComm.wait(cmdTerminateSlave);
                isWaitingCmdTerminateSlave = false;
                
                isTerminate = true;
            }

            // phase3
            JobListType::iterator itEnd = jobList.end();
            for (JobListType::iterator it = jobList.begin(); it != itEnd; ) {
                // 計算
                this->ericalcDH(*it, deltaRho, &deltaHTmp);

                if (isWaitingCmdRequestJob == true) {
                    rComm.wait(cmdRequestJob);
                    isWaitingCmdRequestJob = false;
                }
                rComm.iSendData(cmdRequestJob, root, TAG_REQUEST_JOB);
                isWaitingCmdRequestJob = true;
                
                it = jobList.erase(it);
            }

            // phase4
            if ((isTerminate == true) &&(jobList.empty() == true)) {
                int terminateOK = 0;
                rComm.sendData(terminateOK, root, TAG_TERMINATE_OK);
                break; // break loop
            }
        }
        
        // 後始末
        if (isWaitingCmdAssignJob == true) {
            rComm.cancel(pCmdAssignJob);
            isWaitingCmdAssignJob = false;
        }
        delete[] pCmdAssignJob;
        pCmdAssignJob = NULL;
        
        if (isWaitingCmdTerminateSlave == true) {
            rComm.cancel(cmdTerminateSlave);
            isWaitingCmdTerminateSlave = false;
        }
        
        if (isWaitingCmdRequestJob == true) {
            rComm.wait(cmdRequestJob); // 送信
            isWaitingCmdRequestJob = false;
        }
    }

    //
//     std::cerr << "DfEri_Parallel::getDeltaHpqA_rev2 check1" << std::endl;
    this->loggerTime(" finalize");
    deltaH.mergeSparseMatrix(deltaHTmp);
//     std::cerr << "DfEri_Parallel::getDeltaHpqA_rev2 check2" << std::endl;

    this->cutoffReport();
    this->loggerTime(" end");
//     std::cerr << "DfEri_Parallel::getDeltaHpqA_rev2 check3" << std::endl;
    
    assert(rComm.checkNonBlockingCommunications() == true);
}


