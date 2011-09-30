#include "DfTaskCtrl_Parallel.h"
#include "TlCommunicate.h"

DfTaskCtrl_Parallel::DfTaskCtrl_Parallel(TlSerializeData* pPdfParam) 
    : DfTaskCtrl(pPdfParam) {
    
    const int input_numOfSessions = (*pPdfParam)["num_of_sessions"].getInt();
    this->numOfSessions_ = (input_numOfSessions > 0) ? input_numOfSessions : 1;
}


DfTaskCtrl_Parallel::~DfTaskCtrl_Parallel()
{
}


void DfTaskCtrl_Parallel::cutoffReport()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    rComm.allReduce_SUM(cutoffAll_E1_);
    rComm.allReduce_SUM(cutoffAlive_E1_);
    rComm.allReduce_SUM(cutoffAll_E2_);
    rComm.allReduce_SUM(cutoffAlive_E2_);
    rComm.allReduce_SUM(cutoffAll_schwarz_);
    rComm.allReduce_SUM(cutoffAlive_schwarz_);

    if (rComm.isMaster() == true) {
        DfTaskCtrl::cutoffReport();
    }
    rComm.barrier();
}


bool DfTaskCtrl_Parallel::getQueue(const TlOrbitalInfoObject& orbitalInfo,
                                   const bool isCutoffDistribution,
                                   const int maxGrainSize,
                                   std::vector<Task2>* pTask,
                                   bool initialize)
{
    bool answer = false;
    if (this->isMasterSlave_ == true) {
        answer = this->getQueue_MS(orbitalInfo, isCutoffDistribution,
                                   maxGrainSize, pTask, initialize);
    } else {
        answer = this->getQueue_DC(orbitalInfo, isCutoffDistribution,
                                   maxGrainSize, pTask, initialize);
    }

    return answer;
}


bool DfTaskCtrl_Parallel::getQueue_DC(const TlOrbitalInfoObject& orbitalInfo,
                                      const bool isCutoffDistribution,
                                      const int maxGrainSize,
                                      std::vector<Task2>* pTask,
                                      bool initialize)
{
    assert(pTask != NULL);
    pTask->clear();

    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int numOfProcs = rComm.getNumOfProcs();
    const int globalMaxGrainSize = maxGrainSize * numOfProcs;

    std::vector<Task2> globalTask;
    bool answer = false;
    answer = DfTaskCtrl::getQueue(orbitalInfo, isCutoffDistribution,
                                  globalMaxGrainSize,
                                  &globalTask, initialize);
    if (answer == true) {
        const std::size_t grainSize = globalTask.size();
        const std::size_t localGrainSize = (grainSize + numOfProcs -1) / numOfProcs;
        
        const int rank = rComm.getRank();
        const std::size_t begin = localGrainSize * rank;
        const std::size_t end = std::min(localGrainSize * (rank +1), grainSize);
        pTask->resize(end - begin);
        std::copy(globalTask.begin() + begin,
                  globalTask.begin() + end,
                  pTask->begin());
    }

    return answer;
}


bool DfTaskCtrl_Parallel::getQueue_MS(const TlOrbitalInfoObject& orbitalInfo,
                                      const bool isCutoffDistribution,
                                      const int maxGrainSize,
                                      std::vector<Task2>* pTaskList,
                                      bool initialize)
{
    bool answer = true;
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        answer = this->getQueue_MS_master(orbitalInfo, isCutoffDistribution,
                                        maxGrainSize, pTaskList, initialize);
    } else {
        answer = this->getQueue_MS_slave(orbitalInfo, isCutoffDistribution,
                                         maxGrainSize, pTaskList, initialize);
    }

    return answer;
}


bool DfTaskCtrl_Parallel::getQueue_MS_master(const TlOrbitalInfoObject& orbitalInfo,
                                             const bool isCutoffDistribution,
                                             const int maxGrainSize,
                                             std::vector<Task2>* pTaskList,
                                             bool initialize)
{
    bool answer = true;
    TlCommunicate& rComm = TlCommunicate::getInstance();

    int rank = 0;
    
    // 初期化
    static int numOfFinishStateSlaves = 0;
    if (initialize == true) {
        numOfFinishStateSlaves = 0;
        
        answer = DfTaskCtrl::getQueue(orbitalInfo, isCutoffDistribution,
                                      maxGrainSize, pTaskList, true);
        return answer;
    }


    // 終了メッセージの処理
    static bool isRecvFinishMsg = false;
    static int finishMsg = 0;
    if (isRecvFinishMsg != true) {
        rComm.iReceiveDataFromAnySource(finishMsg, DFTC_FINISH);
        isRecvFinishMsg = true;
    }
    if ((isRecvFinishMsg == true) &&
        (rComm.test(finishMsg, &rank) == true)) {
        rComm.wait(finishMsg);
        ++numOfFinishStateSlaves;
        isRecvFinishMsg = false;
    }
        
    // リクエストの登録
    static RequestListType requestList;
    static bool isRecvSessionID = false;
    static int sessionID = 0;
    static const int maxContinuousSessionRequests = 5; // 調整用変数
    for (int i = 0; i < maxContinuousSessionRequests; ++i) {
        if (isRecvSessionID != true) {
            rComm.iReceiveDataFromAnySource(sessionID, DFTC_SESSION_ID);
            isRecvSessionID = true;
        }

        if ((isRecvSessionID == true) &&
            (rComm.test(sessionID, &rank) == true)) {
            rComm.wait(sessionID);
            isRecvSessionID = false;
            // std::cerr << TlUtils::format("session accept: rank=%d, session=%d", rank, sessionID)
            //           << std::endl;
            
            Request request;
            request.rank = rank;
            request.sessionID = sessionID;
            request.state.reset();
            requestList.push_back(request);
        } else {
            break;
        }
    }

    // 仕事の割り当て
    RequestListType::iterator itEnd = requestList.end();
    for (RequestListType::iterator it = requestList.begin(); it != itEnd; ++it) {
        if (it->state[Request::TASK_ASSIGNED] != true) {
            std::vector<Task2> tmpTaskList;
            const bool hasTask = DfTaskCtrl::getQueue(orbitalInfo, isCutoffDistribution,
                                                      maxGrainSize, &tmpTaskList);
            if (hasTask == true) {
                const std::size_t numOfTasks = tmpTaskList.size();
                const std::size_t unpackTaskListSize = numOfTasks * 2;
                std::vector<index_type> unpackTaskList(unpackTaskListSize);
#pragma omp parallel for
                for (std::size_t i = 0; i < numOfTasks; ++i) {
                    const Task2& task = tmpTaskList[i];
                    unpackTaskList[i*2   ] = task.shellIndex1;
                    unpackTaskList[i*2 +1] = task.shellIndex2;
                }
                it->taskSize = unpackTaskList.size();
                it->unpackTaskList = unpackTaskList;
                it->state[Request::HAS_TASK] = true;
            } else {
                it->taskSize = 0;
                it->state[Request::HAS_TASK] = false;
            }

            it->state[Request::TASK_ASSIGNED] = true;
        }
    }
    
    // 送信
    for (RequestListType::iterator it = requestList.begin(); it != itEnd; ++it) {
        if (it->state[Request::TASK_ASSIGNED] != true) {
            continue;
        }

        if (it->state[Request::SEND_ANSWER] != true) {
            const int rank = it->rank;
            const int sessionID = it->sessionID;

            if (it->state[Request::HAS_TASK] == true) {
                rComm.iSendData(it->taskSize, rank, DFTC_TASK_SIZE + sessionID);
                rComm.iSendDataX((index_type*)&(it->unpackTaskList[0]),
                                 it->taskSize, rank,
                                 DFTC_TASK + sessionID);
                // std::cerr << TlUtils::format("send_task: rank=%d, session=%d", rank, sessionID)
                //           << std::endl;
            } else {
                it->finish = 9999;
                rComm.iSendData(it->finish, rank, DFTC_FINISH + sessionID);
                // std::cerr << TlUtils::format("send_finish: rank=%d, session=%d", rank, sessionID)
                //           << std::endl;
            }

            it->state[Request::SEND_ANSWER] = true;
        } else if (it->state[Request::SESSION_END] != true) {
            if (it->state[Request::HAS_TASK] == true) {
                if ((it->state[Request::WAIT_TASK_SIZE] != true) &&
                    (rComm.test(it->taskSize) == true)) {
                    rComm.wait(it->taskSize);
                    it->state[Request::WAIT_TASK_SIZE] = true;
                } 
                if ((it->state[Request::WAIT_TASK] != true) &&
                    (rComm.test(&(it->unpackTaskList[0])) == true)) {
                    rComm.wait(&(it->unpackTaskList[0]));
                    it->state[Request::WAIT_TASK] = true;
                }
                
                if ((it->state[Request::WAIT_TASK_SIZE] == true) &&
                    (it->state[Request::WAIT_TASK] == true)) {
                    it->state[Request::SESSION_END] = true;
                }
            } else {
                if (rComm.test(it->finish) == true) {
                    rComm.wait(it->finish);
                    it->state[Request::SESSION_END] = true;
                } 
            }
        }
    }

    // 終了フラグのリクエストを削除する
    for (RequestListType::iterator it = requestList.begin(); it != itEnd; ) {
        if (it->state[Request::SESSION_END] == true) {
            it = requestList.erase(it);
            continue;
        }
        ++it;
    }


    // 終了条件のチェック
    const int numOfSlaves = rComm.getNumOfProcs() - 1;
    if (numOfFinishStateSlaves >= numOfSlaves) {
        answer = false;

        if (isRecvFinishMsg == true) {
            rComm.cancel(&finishMsg);
            isRecvFinishMsg = false;
        }
        if (isRecvSessionID == true) {
            rComm.cancel(&sessionID);
            isRecvSessionID = false;
        }
        
        assert(requestList.empty() == true);
    }    

    return answer;
}


bool DfTaskCtrl_Parallel::getQueue_MS_slave(const TlOrbitalInfoObject& orbitalInfo,
                                            const bool isCutoffDistribution,
                                            const int maxGrainSize,
                                            std::vector<Task2>* pTaskList,
                                            bool initialize)
{
    bool answer = true;
    pTaskList->clear();
    TlCommunicate& rComm = TlCommunicate::getInstance();

    // 初期化
    const int numOfSessions = this->numOfSessions_;
    static bool isAllSessionFinished = false;
    static std::vector<Session> sessions(numOfSessions);
    if (initialize == true) {
        isAllSessionFinished = false;
        for (int i = 0; i < numOfSessions; ++i) {
            sessions[i].state.reset();
        }
        
        answer = DfTaskCtrl::getQueue(orbitalInfo, isCutoffDistribution,
                                      maxGrainSize, pTaskList, true);
        return answer;
    }

    if (isAllSessionFinished == true) {
        return false;
    }
    
    const int master = 0;

    for (int i = 0; i < numOfSessions; ++i) {
        Session& session = sessions[i];
        if (session.state[Session::SESSION_END] == true) {
            continue;
        }

        // 終了メッセージの受信
        if (session.state[Session::RECV_FINISH] != true) {
            rComm.iReceiveData(session.recvFinishMsg, master, DFTC_FINISH + i);
            session.state[Session::RECV_FINISH] = true;
        }
        if ((session.state[Session::RECV_FINISH] == true) &&
            (session.state[Session::SESSION_END] != true)) {
            if (rComm.test(session.recvFinishMsg) == true) {
                rComm.wait(session.recvFinishMsg);
                // session.sendExit = 999;
                // rComm.iSendData(session.sendExit, master, DFTC_FINISH);
                session.state[Session::SESSION_END] = true;
            }
        }
        // if (session.state[Session::WAIT_FINISH] == true) {
        //     if (rComm.test(session.sendExit) == true) {
        //         rComm.wait(session.sendExit);
        //         session.state[Session::SESSION_END] = true;
        //         break;
        //     }
        // }

        // task メッセージの受信
        if (session.state[Session::SEND_SESSION] != true) {
            session.id = i;
            rComm.iSendData(session.id, master, DFTC_SESSION_ID);
            // std::cerr << TlUtils::format("send_session")
            //           << std::endl;
            session.state[Session::SEND_SESSION] = true;
        } else if (session.state[Session::WAIT_SESSION] != true) {
            // std::cerr << TlUtils::format("test_session")
            //           << std::endl;
            if (rComm.test(session.id) == true) {
                rComm.wait(session.id);
                sessions[i].state[Session::WAIT_SESSION] = true;
                // std::cerr << TlUtils::format("wait_session")
                //           << std::endl;
            }
        }

        // task sizeの受信
        if (session.state[Session::RECV_SIZE] != true) {
            rComm.iReceiveData(session.size, master, DFTC_TASK_SIZE + i);
            session.state[Session::RECV_SIZE] = true;
        } else if (session.state[Session::WAIT_SIZE] != true) {
            if (rComm.test(session.size) == true) {
                rComm.wait(session.size);
                session.unpackTaskList.resize(session.size);
                session.state[Session::WAIT_SIZE] = true;
            }
        }

        // task listの受信
        if ((session.state[Session::WAIT_SIZE] == true) &&
            (session.state[Session::RECV_TASK_LIST] != true)) {
            session.unpackTaskList.resize(session.size);
            rComm.iReceiveDataX(&(session.unpackTaskList[0]),
                                session.size, master, DFTC_TASK + i);
            session.state[Session::RECV_TASK_LIST] = true;
        }
        if ((session.state[Session::RECV_TASK_LIST] == true) &&
            (session.state[Session::WAIT_TASK_LIST] != true)) {
            if (rComm.test(&(session.unpackTaskList[0])) == true) {
                rComm.wait(&(session.unpackTaskList[0]));
                session.state[Session::WAIT_TASK_LIST] = true;
                
                const std::size_t taskListSize = session.size / 2;
                pTaskList->resize(taskListSize);
                Task2 task;
#pragma omp parallel for private(task)
                for (std::size_t i = 0; i < taskListSize; ++i) {
                    task.shellIndex1 = session.unpackTaskList[i*2   ];
                    task.shellIndex2 = session.unpackTaskList[i*2 +1];
                    (*pTaskList)[i] = task;
                }

                // task 完了
                session.state[Session::SEND_SESSION] = false;
                session.state[Session::WAIT_SESSION] = false;
                session.state[Session::RECV_SIZE] = false;
                session.state[Session::WAIT_SIZE] = false;
                session.state[Session::RECV_TASK_LIST] = false;
                session.state[Session::WAIT_TASK_LIST] = false;

                return true;
            }
        }
    }

    // queueの終了条件
    {
        bool allFinished = true;
        for (int i = 0; i < numOfSessions; ++i) {
            if (sessions[i].state[Session::SESSION_END] != true) {
                allFinished = false;
                break;
            }
        }
        if (allFinished == true) {
            // std::cerr << TlUtils::format("all finished")
            //           << std::endl;
            const int sendExitMsg = 9999;
            rComm.sendData(sendExitMsg, master, DFTC_FINISH);
            
            for (int i = 0; i < numOfSessions; ++i) {
                if ((sessions[i].state[Session::SEND_SESSION] == true) &&
                    (sessions[i].state[Session::WAIT_SESSION] != true)) {
                    // std::cerr << TlUtils::format("cancel:send_session")
                    //           << std::endl;
                    rComm.cancel(sessions[i].id);
                }
                if ((sessions[i].state[Session::RECV_SIZE] == true) &&
                    (sessions[i].state[Session::WAIT_SIZE] != true)) {
                    // std::cerr << TlUtils::format("cancel:recv_size")
                    //           << std::endl;
                    rComm.cancel(sessions[i].size);
                }
            }

            isAllSessionFinished = true;
            answer = false;
        }
    }
    
    return answer;

//     bool answer = true;
//     pTaskList->clear();
//     TlCommunicate& rComm = TlCommunicate::getInstance();

//     // 初期化
//     const int numOfSessions = 2;
//     static bool isAllSessionFinished = false;
//     static std::vector<Session> sessions(numOfSessions);
//     if (initialize == true) {
//         isAllSessionFinished = false;
//         for (int i = 0; i < numOfSessions; ++i) {
//             sessions[i].state.reset();
//         }
        
//         answer = DfTaskCtrl::getQueue(orbitalInfo, isCutoffDistribution,
//                                       maxGrainSize, pTaskList, true);
//         return answer;
//     }

//     if (isAllSessionFinished == true) {
//         return false;
//     }
    
//     const int master = 0;

//     for (int i = 0; i < numOfSessions; ++i) {
//         Session& session = sessions[i];
//         // std::cerr << TlUtils::format("send_loop: rank=%d, session=%d", rComm.getRank(), i)
//         //           << std::endl;
//         if (session.state[Session::SESSION_END] == true) {
//             continue;
//         }

//         // 終了メッセージの受信
//         if (session.state[Session::RECV_FINISH] != true) {
//             rComm.iReceiveData(session.recvFinishMsg, master, DFTC_FINISH + i);
//             session.state[Session::RECV_FINISH] = true;
//         }
//         if ((session.state[Session::RECV_FINISH] == true) &&
//             (session.state[Session::SESSION_END] != true)) {
//             if (rComm.test(session.recvFinishMsg) == true) {
//                 rComm.wait(session.recvFinishMsg);
//                 // session.sendExit = 999;
//                 // rComm.iSendData(session.sendExit, master, DFTC_FINISH);
//                 session.state[Session::SESSION_END] = true;
//                 std::cerr << TlUtils::format("session end: rank=%d, session=%d", rComm.getRank(), session.id)
//                           << std::endl;
//             }
//         }
//         // if (session.state[Session::WAIT_FINISH] == true) {
//         //     if (rComm.test(session.sendExit) == true) {
//         //         rComm.wait(session.sendExit);
//         //         session.state[Session::SESSION_END] = true;
//         //         break;
//         //     }
//         // }

//         // task メッセージの受信
//         if (session.state[Session::SEND_SESSION] != true) {
//             session.id = i;
//             rComm.iSendData(session.id, master, DFTC_SESSION_ID);
//             std::cerr << TlUtils::format("send_session: rank=%d, session=%d", rComm.getRank(), session.id)
//                       << std::endl;
//             session.state[Session::SEND_SESSION] = true;
//         }
//         if ((session.state[Session::SEND_SESSION] == true) &&
//             (session.state[Session::WAIT_SESSION] != true)) {
//             // std::cerr << TlUtils::format("test_session")
//             //           << std::endl;
//             if (rComm.test(session.id) == true) {
//                 rComm.wait(session.id);
//                 sessions[i].state[Session::WAIT_SESSION] = true;
//                 // std::cerr << TlUtils::format("wait_session")
//                 //           << std::endl;
//             }
//         }

//         // task sizeの受信
//         if (session.state[Session::RECV_SIZE] != true) {
//             rComm.iReceiveData(session.size, master, DFTC_TASK_SIZE + i);
//             session.state[Session::RECV_SIZE] = true;
//         } else if (session.state[Session::WAIT_SIZE] != true) {
//             if (rComm.test(session.size) == true) {
//                 rComm.wait(session.size);
//                 session.unpackTaskList.resize(session.size);
//                 session.state[Session::WAIT_SIZE] = true;
//                 std::cerr << TlUtils::format("recv_tasksize: rank=%d, session=%d", rComm.getRank(), session.id)
//                           << std::endl;
//             }
//         }

//         // task listの受信
//         if ((session.state[Session::WAIT_SIZE] == true) &&
//             (session.state[Session::RECV_TASK_LIST] != true)) {
//             session.unpackTaskList.resize(session.size);
//             rComm.iReceiveDataX(&(session.unpackTaskList[0]),
//                                 session.size, master, DFTC_TASK + i);
//             session.state[Session::RECV_TASK_LIST] = true;
//         }
//         if ((session.state[Session::RECV_TASK_LIST] == true) &&
//             (session.state[Session::WAIT_TASK_LIST] != true)) {
//             if (rComm.test(&(session.unpackTaskList[0])) == true) {
//                 rComm.wait(&(session.unpackTaskList[0]));
//                 session.state[Session::WAIT_TASK_LIST] = true;
                
//                 const std::size_t taskListSize = session.size / 2;
//                 pTaskList->resize(taskListSize);
//                 Task2 task;
// #pragma omp parallel for private(task)
//                 for (std::size_t i = 0; i < taskListSize; ++i) {
//                     task.shellIndex1 = session.unpackTaskList[i*2   ];
//                     task.shellIndex2 = session.unpackTaskList[i*2 +1];
//                     (*pTaskList)[i] = task;
//                 }

//                 // task 完了
//                 session.state[Session::SEND_SESSION] = false;
//                 session.state[Session::WAIT_SESSION] = false;
//                 session.state[Session::RECV_SIZE] = false;
//                 session.state[Session::WAIT_SIZE] = false;
//                 session.state[Session::RECV_TASK_LIST] = false;
//                 session.state[Session::WAIT_TASK_LIST] = false;
//                 std::cerr << TlUtils::format("task: rank=%d, session=%d", rComm.getRank(), session.id)
//                           << std::endl;
//             }
//         }
//     }

//     // queueの終了条件
//     {
//         bool allFinished = true;
//         for (int i = 0; i < numOfSessions; ++i) {
//             if (sessions[i].state[Session::SESSION_END] != true) {
//                 allFinished = false;
//                 break;
//             }
//         }
//         if (allFinished == true) {
//             std::cerr << TlUtils::format("all finished")
//                       << std::endl;
//             const int sendExitMsg = 9999;
//             rComm.sendData(sendExitMsg, master, DFTC_FINISH);
            
//             for (int i = 0; i < numOfSessions; ++i) {
//                 if ((sessions[i].state[Session::SEND_SESSION] == true) &&
//                     (sessions[i].state[Session::WAIT_SESSION] != true)) {
//                     // std::cerr << TlUtils::format("cancel:send_session")
//                     //           << std::endl;
//                     rComm.cancel(sessions[i].id);
//                 }
//                 if ((sessions[i].state[Session::RECV_SIZE] == true) &&
//                     (sessions[i].state[Session::WAIT_SIZE] != true)) {
//                     // std::cerr << TlUtils::format("cancel:recv_size")
//                     //           << std::endl;
//                     rComm.cancel(sessions[i].size);
//                 }
//             }

//             isAllSessionFinished = true;
//             answer = false;
//         }
//     }
    
//     return answer;
}


bool DfTaskCtrl_Parallel::getQueue4(const TlOrbitalInfoObject& orbitalInfo,
                                    const TlSparseSymmetricMatrix& schwarzTable,
                                    const int maxGrainSize,
                                    std::vector<Task4>* pTaskList,
                                    bool initialize)
{
    bool answer = false;
    if (this->isMasterSlave_ == true) {
        answer = this->getQueue4_MS(orbitalInfo,
                                    schwarzTable,
                                    maxGrainSize, pTaskList, initialize);
    } else {
        answer = this->getQueue4_DC(orbitalInfo,
                                    schwarzTable,
                                    maxGrainSize, pTaskList, initialize);
    }

    return answer;
}


bool DfTaskCtrl_Parallel::getQueue4_DC(const TlOrbitalInfoObject& orbitalInfo,
                                       const TlSparseSymmetricMatrix& schwarzTable,
                                       const int maxGrainSize,
                                       std::vector<Task4>* pTaskList,
                                       bool initialize)
{
    assert(pTaskList != NULL);
    pTaskList->clear();

    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int globalMaxGrainSize = maxGrainSize * rComm.getNumOfProcs();

    std::vector<Task4> globalTask;
    bool answer = false;
    answer = DfTaskCtrl::getQueue4(orbitalInfo,
                                   schwarzTable,
                                   globalMaxGrainSize,
                                   &globalTask, initialize);
    if (answer == true) {
        const int rank = rComm.getRank();
        std::vector<Task4>::const_iterator begin = std::min(globalTask.begin()
                                                            + maxGrainSize * rank,
                                                            globalTask.end());
        std::vector<Task4>::const_iterator end = std::min(globalTask.begin()
                                                          + maxGrainSize * (rank +1),
                                                          globalTask.end());
        const std::size_t size = std::distance(begin, end);
        pTaskList->resize(size);
        std::copy(begin, end, pTaskList->begin());
    }

    return answer;
}


bool DfTaskCtrl_Parallel::getQueue4_MS(const TlOrbitalInfoObject& orbitalInfo,
                                       const TlSparseSymmetricMatrix& schwarzTable,
                                       const int maxGrainSize,
                                       std::vector<Task4>* pTaskList,
                                       bool initialize)
{
    bool answer = true;
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        answer = this->getQueue4_MS_master(orbitalInfo, schwarzTable,
                                           maxGrainSize, pTaskList, initialize);
    } else {
        answer = this->getQueue4_MS_slave(orbitalInfo, schwarzTable,
                                          maxGrainSize, pTaskList, initialize);
    }

    return answer;
}


bool DfTaskCtrl_Parallel::getQueue4_MS_master(const TlOrbitalInfoObject& orbitalInfo,
                                              const TlSparseSymmetricMatrix& schwarzTable,
                                              const int maxGrainSize,
                                              std::vector<Task4>* pTaskList,
                                              bool initialize)
{
    bool answer = true;
    TlCommunicate& rComm = TlCommunicate::getInstance();

    int rank = 0;
    
    // 初期化
    static int numOfFinishStateSlaves = 0;
    if (initialize == true) {
        numOfFinishStateSlaves = 0;
        
        answer = DfTaskCtrl::getQueue4(orbitalInfo, schwarzTable,
                                       maxGrainSize, pTaskList, true);
        return answer;
    }


    // 終了メッセージの処理
    static bool isRecvFinishMsg = false;
    static int finishMsg = 0;
    if (isRecvFinishMsg != true) {
        rComm.iReceiveDataFromAnySource(finishMsg, DFTC_FINISH);
        isRecvFinishMsg = true;
    }
    if ((isRecvFinishMsg == true) &&
        (rComm.test(finishMsg, &rank) == true)) {
        rComm.wait(finishMsg);
        ++numOfFinishStateSlaves;
        isRecvFinishMsg = false;
    }
        
    // リクエストの登録
    static RequestListType requestList;
    static bool isRecvSessionID = false;
    static int sessionID = 0;
    static const int maxContinuousSessionRequests = 5; // 調整用変数
    for (int i = 0; i < maxContinuousSessionRequests; ++i) {
        if (isRecvSessionID != true) {
            rComm.iReceiveDataFromAnySource(sessionID, DFTC_SESSION_ID);
            isRecvSessionID = true;
        }

        if ((isRecvSessionID == true) &&
            (rComm.test(sessionID, &rank) == true)) {
            rComm.wait(sessionID);
            isRecvSessionID = false;

            Request request;
            request.rank = rank;
            request.sessionID = sessionID;
            request.state.reset();
            requestList.push_back(request);
        } else {
            break;
        }
    }

    // 仕事の割り当て
    RequestListType::iterator itEnd = requestList.end();
    for (RequestListType::iterator it = requestList.begin(); it != itEnd; ++it) {
        if (it->state[Request::TASK_ASSIGNED] != true) {
            std::vector<Task4> tmpTaskList;
            const bool hasTask = DfTaskCtrl::getQueue4(orbitalInfo, schwarzTable,
                                                       maxGrainSize, &tmpTaskList);
            if (hasTask == true) {
                const std::size_t numOfTasks = tmpTaskList.size();
                const std::size_t unpackTaskListSize = numOfTasks * 4;
                std::vector<index_type> unpackTaskList(unpackTaskListSize);
#pragma omp parallel for
                for (std::size_t i = 0; i < numOfTasks; ++i) {
                    const Task4& task = tmpTaskList[i];
                    unpackTaskList[i*4   ] = task.shellIndex1;
                    unpackTaskList[i*4 +1] = task.shellIndex2;
                    unpackTaskList[i*4 +2] = task.shellIndex3;
                    unpackTaskList[i*4 +3] = task.shellIndex4;
                }
                it->taskSize = unpackTaskList.size();
                it->unpackTaskList = unpackTaskList;
                it->state[Request::HAS_TASK] = true;
            } else {
                it->taskSize = 0;
                it->state[Request::HAS_TASK] = false;
            }

            it->state[Request::TASK_ASSIGNED] = true;
        }
    }
    
    // 送信
    for (RequestListType::iterator it = requestList.begin(); it != itEnd; ++it) {
        if (it->state[Request::TASK_ASSIGNED] != true) {
            continue;
        }

        if (it->state[Request::SEND_ANSWER] != true) {
            const int rank = it->rank;
            const int sessionID = it->sessionID;

            if (it->state[Request::HAS_TASK] == true) {
                rComm.iSendData(it->taskSize, rank, DFTC_TASK_SIZE + sessionID);
                rComm.iSendDataX((index_type*)&(it->unpackTaskList[0]), it->taskSize, rank, DFTC_TASK + sessionID);
            } else {
                it->finish = 9999;
                rComm.iSendData(it->finish, rank, DFTC_FINISH + sessionID);
            }

            it->state[Request::SEND_ANSWER] = true;
        } else if (it->state[Request::SESSION_END] != true) {
            if (it->state[Request::HAS_TASK] == true) {
                if ((it->state[Request::WAIT_TASK_SIZE] != true) &&
                    (rComm.test(it->taskSize) == true)) {
                    rComm.wait(it->taskSize);
                    it->state[Request::WAIT_TASK_SIZE] = true;
                } 
                if ((it->state[Request::WAIT_TASK] != true) &&
                    (rComm.test(&(it->unpackTaskList[0])) == true)) {
                    rComm.wait(&(it->unpackTaskList[0]));
                    it->state[Request::WAIT_TASK] = true;
                }
                
                if ((it->state[Request::WAIT_TASK_SIZE] == true) &&
                    (it->state[Request::WAIT_TASK] == true)) {
                    it->state[Request::SESSION_END] = true;
                }
            } else {
                if (rComm.test(it->finish) == true) {
                    rComm.wait(it->finish);
                    it->state[Request::SESSION_END] = true;
                } 
            }
        }
    }

    // 終了フラグのリクエストを削除する
    for (RequestListType::iterator it = requestList.begin(); it != itEnd; ) {
        if (it->state[Request::SESSION_END] == true) {
            it = requestList.erase(it);
            continue;
        }
        ++it;
    }


    // 終了条件のチェック
    const int numOfSlaves = rComm.getNumOfProcs() - 1;
    if (numOfFinishStateSlaves >= numOfSlaves) {
        answer = false;

        if (isRecvFinishMsg == true) {
            rComm.cancel(&finishMsg);
            isRecvFinishMsg = false;
        }
        if (isRecvSessionID == true) {
            rComm.cancel(&sessionID);
            isRecvSessionID = false;
        }
        
        assert(requestList.empty() == true);
    }    

    return answer;
}


bool DfTaskCtrl_Parallel::getQueue4_MS_slave(const TlOrbitalInfoObject& orbitalInfo,
                                             const TlSparseSymmetricMatrix& schwarzTable,
                                             const int maxGrainSize,
                                             std::vector<Task4>* pTaskList,
                                             bool initialize)
{
    bool answer = true;
    pTaskList->clear();
    TlCommunicate& rComm = TlCommunicate::getInstance();

    // 初期化
    const int numOfSessions = this->numOfSessions_;
    static bool isAllSessionFinished = false;
    static std::vector<Session> sessions(numOfSessions);
    if (initialize == true) {
        isAllSessionFinished = false;
        for (int i = 0; i < numOfSessions; ++i) {
            sessions[i].state.reset();
        }
        
        answer = DfTaskCtrl::getQueue4(orbitalInfo, schwarzTable,
                                       maxGrainSize, pTaskList, true);
        return answer;
    }

    if (isAllSessionFinished == true) {
        return false;
    }
    
    const int master = 0;

    for (int i = 0; i < numOfSessions; ++i) {
        Session& session = sessions[i];
        if (session.state[Session::SESSION_END] == true) {
            continue;
        }

        // 終了メッセージの受信
        if (session.state[Session::RECV_FINISH] != true) {
            rComm.iReceiveData(session.recvFinishMsg, master, DFTC_FINISH + i);
            session.state[Session::RECV_FINISH] = true;
        }
        if ((session.state[Session::RECV_FINISH] == true) &&
            (session.state[Session::SESSION_END] != true)) {
            if (rComm.test(session.recvFinishMsg) == true) {
                rComm.wait(session.recvFinishMsg);
                // session.sendExit = 999;
                // rComm.iSendData(session.sendExit, master, DFTC_FINISH);
                session.state[Session::SESSION_END] = true;
            }
        }
        // if (session.state[Session::WAIT_FINISH] == true) {
        //     if (rComm.test(session.sendExit) == true) {
        //         rComm.wait(session.sendExit);
        //         session.state[Session::SESSION_END] = true;
        //         break;
        //     }
        // }

        // task メッセージの受信
        if (session.state[Session::SEND_SESSION] != true) {
            session.id = i;
            rComm.iSendData(session.id, master, DFTC_SESSION_ID);
            // std::cerr << TlUtils::format("send_session")
            //           << std::endl;
            session.state[Session::SEND_SESSION] = true;
        } else if (session.state[Session::WAIT_SESSION] != true) {
            // std::cerr << TlUtils::format("test_session")
            //           << std::endl;
            if (rComm.test(session.id) == true) {
                rComm.wait(session.id);
                sessions[i].state[Session::WAIT_SESSION] = true;
                // std::cerr << TlUtils::format("wait_session")
                //           << std::endl;
            }
        }

        // task sizeの受信
        if (session.state[Session::RECV_SIZE] != true) {
            rComm.iReceiveData(session.size, master, DFTC_TASK_SIZE + i);
            session.state[Session::RECV_SIZE] = true;
        } else if (session.state[Session::WAIT_SIZE] != true) {
            if (rComm.test(session.size) == true) {
                rComm.wait(session.size);
                session.unpackTaskList.resize(session.size);
                session.state[Session::WAIT_SIZE] = true;
            }
        }

        // task listの受信
        if ((session.state[Session::WAIT_SIZE] == true) &&
            (session.state[Session::RECV_TASK_LIST] != true)) {
            session.unpackTaskList.resize(session.size);
            rComm.iReceiveDataX(&(session.unpackTaskList[0]),
                                session.size, master, DFTC_TASK + i);
            session.state[Session::RECV_TASK_LIST] = true;
        }
        if ((session.state[Session::RECV_TASK_LIST] == true) &&
            (session.state[Session::WAIT_TASK_LIST] != true)) {
            if (rComm.test(&(session.unpackTaskList[0])) == true) {
                rComm.wait(&(session.unpackTaskList[0]));
                session.state[Session::WAIT_TASK_LIST] = true;
                
                const std::size_t taskListSize = session.size / 4;
                pTaskList->resize(taskListSize);
                Task4 task;
#pragma omp parallel for private(task)
                for (std::size_t i = 0; i < taskListSize; ++i) {
                    task.shellIndex1 = session.unpackTaskList[i*4   ];
                    task.shellIndex2 = session.unpackTaskList[i*4 +1];
                    task.shellIndex3 = session.unpackTaskList[i*4 +2];
                    task.shellIndex4 = session.unpackTaskList[i*4 +3];
                    (*pTaskList)[i] = task;
                }

                // task 完了
                session.state[Session::SEND_SESSION] = false;
                session.state[Session::WAIT_SESSION] = false;
                session.state[Session::RECV_SIZE] = false;
                session.state[Session::WAIT_SIZE] = false;
                session.state[Session::RECV_TASK_LIST] = false;
                session.state[Session::WAIT_TASK_LIST] = false;

                return true;
            }
        }
    }

    // queueの終了条件
    {
        bool allFinished = true;
        for (int i = 0; i < numOfSessions; ++i) {
            if (sessions[i].state[Session::SESSION_END] != true) {
                allFinished = false;
                break;
            }
        }
        if (allFinished == true) {
            // std::cerr << TlUtils::format("all finished")
            //           << std::endl;
            const int sendExitMsg = 9999;
            rComm.sendData(sendExitMsg, master, DFTC_FINISH);
            
            for (int i = 0; i < numOfSessions; ++i) {
                if ((sessions[i].state[Session::SEND_SESSION] == true) &&
                    (sessions[i].state[Session::WAIT_SESSION] != true)) {
                    // std::cerr << TlUtils::format("cancel:send_session")
                    //           << std::endl;
                    rComm.cancel(sessions[i].id);
                }
                if ((sessions[i].state[Session::RECV_SIZE] == true) &&
                    (sessions[i].state[Session::WAIT_SIZE] != true)) {
                    // std::cerr << TlUtils::format("cancel:recv_size")
                    //           << std::endl;
                    rComm.cancel(sessions[i].size);
                }
            }

            isAllSessionFinished = true;
            answer = false;
        }
    }
    
    return answer;
}


bool DfTaskCtrl_Parallel::getQueue_Force4(const TlOrbitalInfoObject& orbitalInfo,
                                          const TlSparseSymmetricMatrix& schwarzTable,
                                          const int maxGrainSize,
                                          std::vector<Task4>* pTaskList,
                                          bool initialize)
{
    bool answer = false;
    if (this->isMasterSlave_ == true) {
        answer = this->getQueue_Force4_MS(orbitalInfo,
                                          schwarzTable,
                                          maxGrainSize,
                                          pTaskList, initialize);
    } else {
        answer = this->getQueue_Force4_DC(orbitalInfo,
                                          schwarzTable,
                                          maxGrainSize,
                                          pTaskList, initialize);
    }

    return answer;
}


bool DfTaskCtrl_Parallel::getQueue_Force4_DC(const TlOrbitalInfoObject& orbitalInfo,
                                             const TlSparseSymmetricMatrix& schwarzTable,
                                             const int maxGrainSize,
                                             std::vector<Task4>* pTaskList,
                                             bool initialize)
{
    assert(pTaskList != NULL);
    pTaskList->clear();

    TlCommunicate& rComm = TlCommunicate::getInstance();
    const int globalMaxGrainSize = maxGrainSize * rComm.getNumOfProcs();

    std::vector<Task4> globalTask;
    bool answer = false;
    answer = DfTaskCtrl::getQueue_Force4(orbitalInfo,
                                         schwarzTable,
                                         globalMaxGrainSize,
                                         &globalTask,
                                         initialize);
    if (answer == true) {
        const int rank = rComm.getRank();
        std::vector<Task4>::const_iterator begin = std::min(globalTask.begin()
                                                            + maxGrainSize * rank,
                                                            globalTask.end());
        std::vector<Task4>::const_iterator end = std::min(globalTask.begin()
                                                          + maxGrainSize * (rank +1),
                                                          globalTask.end());
        const std::size_t size = std::distance(begin, end);
        pTaskList->resize(size);
        std::copy(begin, end, pTaskList->begin());
    }

    return answer;
}


bool DfTaskCtrl_Parallel::getQueue_Force4_MS(const TlOrbitalInfoObject& orbitalInfo,
                                             const TlSparseSymmetricMatrix& schwarzTable,
                                             const int maxGrainSize,
                                             std::vector<Task4>* pTaskList,
                                             bool initialize)
{
    bool answer = true;
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        answer = this->getQueue_Force4_MS_master(orbitalInfo, schwarzTable,
                                                 maxGrainSize, pTaskList, initialize);
    } else {
        answer = this->getQueue_Force4_MS_slave(orbitalInfo, schwarzTable,
                                                maxGrainSize, pTaskList, initialize);
    }

    return answer;
}


bool DfTaskCtrl_Parallel::getQueue_Force4_MS_master(const TlOrbitalInfoObject& orbitalInfo,
                                                    const TlSparseSymmetricMatrix& schwarzTable,
                                                    const int maxGrainSize,
                                                    std::vector<Task4>* pTaskList,
                                                    bool initialize)
{
    bool answer = true;
    TlCommunicate& rComm = TlCommunicate::getInstance();

    int rank = 0;
    
    // 初期化
    static int numOfFinishStateSlaves = 0;
    if (initialize == true) {
        numOfFinishStateSlaves = 0;
        
        answer = DfTaskCtrl::getQueue_Force4(orbitalInfo, schwarzTable,
                                             maxGrainSize, pTaskList, true);
        return answer;
    }


    // 終了メッセージの処理
    static bool isRecvFinishMsg = false;
    static int finishMsg = 0;
    if (isRecvFinishMsg != true) {
        rComm.iReceiveDataFromAnySource(finishMsg, DFTC_FINISH);
        isRecvFinishMsg = true;
    }
    if ((isRecvFinishMsg == true) &&
        (rComm.test(finishMsg, &rank) == true)) {
        rComm.wait(finishMsg);
        ++numOfFinishStateSlaves;
        isRecvFinishMsg = false;
    }
        
    // リクエストの登録
    static RequestListType requestList;
    static bool isRecvSessionID = false;
    static int sessionID = 0;
    static const int maxContinuousSessionRequests = 5; // 調整用変数
    for (int i = 0; i < maxContinuousSessionRequests; ++i) {
        if (isRecvSessionID != true) {
            rComm.iReceiveDataFromAnySource(sessionID, DFTC_SESSION_ID);
            isRecvSessionID = true;
        }

        if ((isRecvSessionID == true) &&
            (rComm.test(sessionID, &rank) == true)) {
            rComm.wait(sessionID);
            isRecvSessionID = false;

            Request request;
            request.rank = rank;
            request.sessionID = sessionID;
            request.state.reset();
            requestList.push_back(request);
        } else {
            break;
        }
    }

    // 仕事の割り当て
    RequestListType::iterator itEnd = requestList.end();
    for (RequestListType::iterator it = requestList.begin(); it != itEnd; ++it) {
        if (it->state[Request::TASK_ASSIGNED] != true) {
            std::vector<Task4> tmpTaskList;
            const bool hasTask = DfTaskCtrl::getQueue_Force4(orbitalInfo, schwarzTable,
                                                             maxGrainSize, &tmpTaskList);
            if (hasTask == true) {
                const std::size_t numOfTasks = tmpTaskList.size();
                const std::size_t unpackTaskListSize = numOfTasks * 4;
                std::vector<index_type> unpackTaskList(unpackTaskListSize);
#pragma omp parallel for
                for (std::size_t i = 0; i < numOfTasks; ++i) {
                    const Task4& task = tmpTaskList[i];
                    unpackTaskList[i*4   ] = task.shellIndex1;
                    unpackTaskList[i*4 +1] = task.shellIndex2;
                    unpackTaskList[i*4 +2] = task.shellIndex3;
                    unpackTaskList[i*4 +3] = task.shellIndex4;
                }
                it->taskSize = unpackTaskList.size();
                it->unpackTaskList = unpackTaskList;
                it->state[Request::HAS_TASK] = true;
            } else {
                it->taskSize = 0;
                it->state[Request::HAS_TASK] = false;
            }

            it->state[Request::TASK_ASSIGNED] = true;
        }
    }
    
    // 送信
    for (RequestListType::iterator it = requestList.begin(); it != itEnd; ++it) {
        if (it->state[Request::TASK_ASSIGNED] != true) {
            continue;
        }

        if (it->state[Request::SEND_ANSWER] != true) {
            const int rank = it->rank;
            const int sessionID = it->sessionID;

            if (it->state[Request::HAS_TASK] == true) {
                rComm.iSendData(it->taskSize, rank, DFTC_TASK_SIZE + sessionID);
                rComm.iSendDataX((index_type*)&(it->unpackTaskList[0]), it->taskSize, rank, DFTC_TASK + sessionID);
            } else {
                it->finish = 9999;
                rComm.iSendData(it->finish, rank, DFTC_FINISH + sessionID);
            }

            it->state[Request::SEND_ANSWER] = true;
        } else if (it->state[Request::SESSION_END] != true) {
            if (it->state[Request::HAS_TASK] == true) {
                if ((it->state[Request::WAIT_TASK_SIZE] != true) &&
                    (rComm.test(it->taskSize) == true)) {
                    rComm.wait(it->taskSize);
                    it->state[Request::WAIT_TASK_SIZE] = true;
                } 
                if ((it->state[Request::WAIT_TASK] != true) &&
                    (rComm.test(&(it->unpackTaskList[0])) == true)) {
                    rComm.wait(&(it->unpackTaskList[0]));
                    it->state[Request::WAIT_TASK] = true;
                }
                
                if ((it->state[Request::WAIT_TASK_SIZE] == true) &&
                    (it->state[Request::WAIT_TASK] == true)) {
                    it->state[Request::SESSION_END] = true;
                }
            } else {
                if (rComm.test(it->finish) == true) {
                    rComm.wait(it->finish);
                    it->state[Request::SESSION_END] = true;
                } 
            }
        }
    }

    // 終了フラグのリクエストを削除する
    for (RequestListType::iterator it = requestList.begin(); it != itEnd; ) {
        if (it->state[Request::SESSION_END] == true) {
            it = requestList.erase(it);
            continue;
        }
        ++it;
    }


    // 終了条件のチェック
    const int numOfSlaves = rComm.getNumOfProcs() - 1;
    if (numOfFinishStateSlaves >= numOfSlaves) {
        answer = false;

        if (isRecvFinishMsg == true) {
            rComm.cancel(&finishMsg);
            isRecvFinishMsg = false;
        }
        if (isRecvSessionID == true) {
            rComm.cancel(&sessionID);
            isRecvSessionID = false;
        }
        
        assert(requestList.empty() == true);
    }    

    return answer;
}


bool DfTaskCtrl_Parallel::getQueue_Force4_MS_slave(const TlOrbitalInfoObject& orbitalInfo,
                                                   const TlSparseSymmetricMatrix& schwarzTable,
                                                   const int maxGrainSize,
                                                   std::vector<Task4>* pTaskList,
                                                   bool initialize)
{
    bool answer = true;
    pTaskList->clear();
    TlCommunicate& rComm = TlCommunicate::getInstance();

    // 初期化
    const int numOfSessions = this->numOfSessions_;
    static bool isAllSessionFinished = false;
    static std::vector<Session> sessions(numOfSessions);
    if (initialize == true) {
        isAllSessionFinished = false;
        for (int i = 0; i < numOfSessions; ++i) {
            sessions[i].state.reset();
        }
        
        answer = DfTaskCtrl::getQueue_Force4(orbitalInfo, schwarzTable,
                                             maxGrainSize, pTaskList, true);
        return answer;
    }

    if (isAllSessionFinished == true) {
        return false;
    }
    
    const int master = 0;

    for (int i = 0; i < numOfSessions; ++i) {
        Session& session = sessions[i];
        if (session.state[Session::SESSION_END] == true) {
            continue;
        }

        // 終了メッセージの受信
        if (session.state[Session::RECV_FINISH] != true) {
            rComm.iReceiveData(session.recvFinishMsg, master, DFTC_FINISH + i);
            session.state[Session::RECV_FINISH] = true;
        }
        if ((session.state[Session::RECV_FINISH] == true) &&
            (session.state[Session::SESSION_END] != true)) {
            if (rComm.test(session.recvFinishMsg) == true) {
                rComm.wait(session.recvFinishMsg);
                // session.sendExit = 999;
                // rComm.iSendData(session.sendExit, master, DFTC_FINISH);
                session.state[Session::SESSION_END] = true;
            }
        }
        // if (session.state[Session::WAIT_FINISH] == true) {
        //     if (rComm.test(session.sendExit) == true) {
        //         rComm.wait(session.sendExit);
        //         session.state[Session::SESSION_END] = true;
        //         break;
        //     }
        // }

        // task メッセージの受信
        if (session.state[Session::SEND_SESSION] != true) {
            session.id = i;
            rComm.iSendData(session.id, master, DFTC_SESSION_ID);
            // std::cerr << TlUtils::format("send_session")
            //           << std::endl;
            session.state[Session::SEND_SESSION] = true;
        } else if (session.state[Session::WAIT_SESSION] != true) {
            // std::cerr << TlUtils::format("test_session")
            //           << std::endl;
            if (rComm.test(session.id) == true) {
                rComm.wait(session.id);
                sessions[i].state[Session::WAIT_SESSION] = true;
                // std::cerr << TlUtils::format("wait_session")
                //           << std::endl;
            }
        }

        // task sizeの受信
        if (session.state[Session::RECV_SIZE] != true) {
            rComm.iReceiveData(session.size, master, DFTC_TASK_SIZE + i);
            session.state[Session::RECV_SIZE] = true;
        } else if (session.state[Session::WAIT_SIZE] != true) {
            if (rComm.test(session.size) == true) {
                rComm.wait(session.size);
                session.unpackTaskList.resize(session.size);
                session.state[Session::WAIT_SIZE] = true;
            }
        }

        // task listの受信
        if ((session.state[Session::WAIT_SIZE] == true) &&
            (session.state[Session::RECV_TASK_LIST] != true)) {
            session.unpackTaskList.resize(session.size);
            rComm.iReceiveDataX(&(session.unpackTaskList[0]),
                                session.size, master, DFTC_TASK + i);
            session.state[Session::RECV_TASK_LIST] = true;
        }
        if ((session.state[Session::RECV_TASK_LIST] == true) &&
            (session.state[Session::WAIT_TASK_LIST] != true)) {
            if (rComm.test(&(session.unpackTaskList[0])) == true) {
                rComm.wait(&(session.unpackTaskList[0]));
                session.state[Session::WAIT_TASK_LIST] = true;
                
                const std::size_t taskListSize = session.size / 4;
                pTaskList->resize(taskListSize);
                Task4 task;
#pragma omp parallel for private(task)
                for (std::size_t i = 0; i < taskListSize; ++i) {
                    task.shellIndex1 = session.unpackTaskList[i*4   ];
                    task.shellIndex2 = session.unpackTaskList[i*4 +1];
                    task.shellIndex3 = session.unpackTaskList[i*4 +2];
                    task.shellIndex4 = session.unpackTaskList[i*4 +3];
                    (*pTaskList)[i] = task;
                }

                // task 完了
                session.state[Session::SEND_SESSION] = false;
                session.state[Session::WAIT_SESSION] = false;
                session.state[Session::RECV_SIZE] = false;
                session.state[Session::WAIT_SIZE] = false;
                session.state[Session::RECV_TASK_LIST] = false;
                session.state[Session::WAIT_TASK_LIST] = false;

                return true;
            }
        }
    }

    // queueの終了条件
    {
        bool allFinished = true;
        for (int i = 0; i < numOfSessions; ++i) {
            if (sessions[i].state[Session::SESSION_END] != true) {
                allFinished = false;
                break;
            }
        }
        if (allFinished == true) {
            // std::cerr << TlUtils::format("all finished")
            //           << std::endl;
            const int sendExitMsg = 9999;
            rComm.sendData(sendExitMsg, master, DFTC_FINISH);
            
            for (int i = 0; i < numOfSessions; ++i) {
                if ((sessions[i].state[Session::SEND_SESSION] == true) &&
                    (sessions[i].state[Session::WAIT_SESSION] != true)) {
                    // std::cerr << TlUtils::format("cancel:send_session")
                    //           << std::endl;
                    rComm.cancel(sessions[i].id);
                }
                if ((sessions[i].state[Session::RECV_SIZE] == true) &&
                    (sessions[i].state[Session::WAIT_SIZE] != true)) {
                    // std::cerr << TlUtils::format("cancel:recv_size")
                    //           << std::endl;
                    rComm.cancel(sessions[i].size);
                }
            }

            isAllSessionFinished = true;
            answer = false;
        }
    }
    
    return answer;
}
