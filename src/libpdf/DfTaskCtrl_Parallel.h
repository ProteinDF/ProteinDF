// Copyright (C) 2002-2014 The ProteinDF project
// see also AUTHORS and README.
//
// This file is part of ProteinDF.
//
// ProteinDF is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ProteinDF is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

#ifndef DFTASKCTRL_PARALLEL_H
#define DFTASKCTRL_PARALLEL_H

#include <bitset>
#include <list>
#include "DfTaskCtrl.h"

class DfTaskCtrl_Parallel : public DfTaskCtrl {
   public:
    DfTaskCtrl_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfTaskCtrl_Parallel();

   public:
    virtual bool getQueue(const std::size_t maxIndeces,
                          const std::size_t maxGrainSize,
                          std::vector<std::size_t>* pTasks,
                          bool initialize = false);

    virtual bool getQueue(const TlOrbitalInfoObject& orbitalInfo,
                          const int maxGrainSize, std::vector<Task>* pTask,
                          bool initialize = false);

    virtual bool getQueue2(const TlOrbitalInfoObject& orbitalInfo,
                           const bool isCutoffByDistibution,
                           const int maxGrainSize, std::vector<Task2>* pTask,
                           bool initialize = false);

    /// 2つのインデックスのタスクを返す(軌道情報が異なる場合)
    virtual bool getQueue2(const TlOrbitalInfoObject& orbitalInfo1,
                           const TlOrbitalInfoObject& orbitalInfo2,
                           const bool isCutoffByDistibution,
                           const int maxGrainSize, std::vector<Task2>* pTask,
                           bool initialize = false);

    virtual bool getQueue4(const TlOrbitalInfoObject& orbitalInfo,
                           const TlSparseSymmetricMatrix& schwarzTable,
                           const int maxGrainSize,
                           std::vector<Task4>* pTaskList,
                           bool initialize = false);

    virtual bool getQueue_Force4(const TlOrbitalInfoObject& orbitalInfo,
                                 const TlSparseSymmetricMatrix& schwarzTable,
                                 const int maxGrainSize,
                                 std::vector<Task4>* pTaskList,
                                 bool initialize = false);

   protected:
    enum {
        TCMSG_MASTER_TO_SLAVE = 1001,
        TCMSG_SLAVE_TO_MASTER = 1002,
        TCMSG_TASK_SIZE = 1003,
        TCMSG_TASK = 1004,
    };

   protected:
    enum {
        DFTC_SESSION_ID = 5000,
        DFTC_TASK_SIZE = 5100,
        DFTC_TASK = 5200,
        DFTC_FINISH = 5900
    };

    struct Request {
        enum {
            TASK_ASSIGNED = 0,
            HAS_TASK,
            SEND_ANSWER,
            WAIT_TASK_SIZE,
            WAIT_TASK,
            SESSION_END,
            NUM_OF_REQUEST_STATES
        };

        Request()
            : rank(0),
              sessionID(0),
              taskSize(0),
              unpackTaskList(std::vector<index_type>()),
              finish(0),
              state(0){};

        int rank;
        int sessionID;
        std::size_t taskSize;
        std::vector<index_type> unpackTaskList;
        int finish;
        std::bitset<NUM_OF_REQUEST_STATES> state;
    };
    typedef std::list<Request> RequestListType;

    struct Session {
        enum {
            SEND_SESSION = 0,
            WAIT_SESSION,
            RECV_SIZE,
            WAIT_SIZE,
            RECV_TASK_LIST,
            WAIT_TASK_LIST,
            RECV_FINISH,
            SESSION_END,
            NUM_OF_SESSION_STATES
        };

        int id;
        std::size_t size;
        std::vector<index_type> unpackTaskList;
        int recvFinishMsg;
        int sendExit;
        std::bitset<NUM_OF_SESSION_STATES> state;
    };

   protected:
    virtual void prescreeningReport();
    virtual void cutoffReport();
    void cutoffReport_MS();

   protected:
    bool getQueue_DC(const std::size_t maxIndeces,
                     const std::size_t maxGrainSize,
                     std::vector<std::size_t>* pTasks, bool initialize = false);

    bool getQueue_DC(const TlOrbitalInfoObject& orbitalInfo,
                     const int maxGrainSize, std::vector<Task>* pTask,
                     bool initialize);

    bool getQueue2_DC(const TlOrbitalInfoObject& orbitalInfo,
                      const bool isCutoffByDistibution, const int maxGrainSize,
                      std::vector<Task2>* pTask, bool initialize = false);

    bool getQueue2_MS(const TlOrbitalInfoObject& orbitalInfo,
                      const bool isCutoffByDistibution, const int maxGrainSize,
                      std::vector<Task2>* pTask, bool initialize = false);

    bool getQueue2_MS_master(const TlOrbitalInfoObject& orbitalInfo,
                             const bool isCutoffDistribution,
                             const int maxGrainSize,
                             std::vector<Task2>* pTaskList, bool initialize);

    bool getQueue2_MS_slave(const TlOrbitalInfoObject& orbitalInfo,
                            const bool isCutoffDistribution,
                            const int maxGrainSize,
                            std::vector<Task2>* pTaskList, bool initialize);

    //
    bool getQueue2_DC(const TlOrbitalInfoObject& orbitalInfo1,
                      const TlOrbitalInfoObject& orbitalInfo2,
                      const bool isCutoffByDistibution, const int maxGrainSize,
                      std::vector<Task2>* pTask, bool initialize = false);

    bool getQueue2_MS(const TlOrbitalInfoObject& orbitalInfo1,
                      const TlOrbitalInfoObject& orbitalInfo2,
                      const bool isCutoffByDistibution, const int maxGrainSize,
                      std::vector<Task2>* pTask, bool initialize = false);

    bool getQueue2_MS_master(const TlOrbitalInfoObject& orbitalInfo1,
                             const TlOrbitalInfoObject& orbitalInfo2,
                             const bool isCutoffDistribution,
                             const int maxGrainSize,
                             std::vector<Task2>* pTaskList, bool initialize);

    bool getQueue2_MS_slave(const TlOrbitalInfoObject& orbitalInfo1,
                            const TlOrbitalInfoObject& orbitalInfo2,
                            const bool isCutoffDistribution,
                            const int maxGrainSize,
                            std::vector<Task2>* pTaskList, bool initialize);
    //

    bool getQueue4_DC(const TlOrbitalInfoObject& orbitalInfo,
                      const TlSparseSymmetricMatrix& schwarzTable,
                      const int maxGrainSize, std::vector<Task4>* pTaskList,
                      bool initialize = false);
    bool getQueue4_MS(const TlOrbitalInfoObject& orbitalInfo,
                      const TlSparseSymmetricMatrix& schwarzTable,
                      const int maxGrainSize, std::vector<Task4>* pTaskList,
                      bool initialize = false);
    bool getQueue4_MS_master(const TlOrbitalInfoObject& orbitalInfo,
                             const TlSparseSymmetricMatrix& schwarzTable,
                             const int maxGrainSize,
                             std::vector<Task4>* pTaskList,
                             bool initialize = false);
    bool getQueue4_MS_slave(const TlOrbitalInfoObject& orbitalInfo,
                            const TlSparseSymmetricMatrix& schwarzTable,
                            const int maxGrainSize,
                            std::vector<Task4>* pTaskList,
                            bool initialize = false);

    bool getQueue_Force4_DC(const TlOrbitalInfoObject& orbitalInfo,
                            const TlSparseSymmetricMatrix& schwarzTable,
                            const int maxGrainSize,
                            std::vector<Task4>* pTaskList,
                            bool initialize = false);
    bool getQueue_Force4_MS(const TlOrbitalInfoObject& orbitalInfo,
                            const TlSparseSymmetricMatrix& schwarzTable,
                            const int maxGrainSize,
                            std::vector<Task4>* pTaskList,
                            bool initialize = false);
    bool getQueue_Force4_MS_master(const TlOrbitalInfoObject& orbitalInfo,
                                   const TlSparseSymmetricMatrix& schwarzTable,
                                   const int maxGrainSize,
                                   std::vector<Task4>* pTaskList,
                                   bool initialize = false);
    bool getQueue_Force4_MS_slave(const TlOrbitalInfoObject& orbitalInfo,
                                  const TlSparseSymmetricMatrix& schwarzTable,
                                  const int maxGrainSize,
                                  std::vector<Task4>* pTaskList,
                                  bool initialize = false);

   protected:
    int numOfSessions_;
};

#endif  // DFTASKCTRL_PARALLEL_H
