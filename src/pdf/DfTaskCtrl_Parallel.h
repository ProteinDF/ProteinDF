#ifndef DFTASKCTRL_PARALLEL_H
#define DFTASKCTRL_PARALLEL_H

#include <list>
#include <bitset>
#include "DfTaskCtrl.h"

class DfTaskCtrl_Parallel : public DfTaskCtrl {
public:
    DfTaskCtrl_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfTaskCtrl_Parallel();

public:
    virtual void cutoffReport();
    virtual bool getQueue(const TlOrbitalInfoObject& orbitalInfo,
                          const bool isCutoffByDistibution,
                          const int maxGrainSize,
                          std::vector<Task2>* pTask,
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
        DFTC_TASK_SIZE  = 5100,
        DFTC_TASK       = 5200,
        DFTC_FINISH     = 5900
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
    bool getQueue_DC(const TlOrbitalInfoObject& orbitalInfo,
                     const bool isCutoffByDistibution,
                     const int maxGrainSize,
                     std::vector<Task2>* pTask,
                     bool initialize = false);

    bool getQueue_MS(const TlOrbitalInfoObject& orbitalInfo,
                     const bool isCutoffByDistibution,
                     const int maxGrainSize,
                     std::vector<Task2>* pTask,
                     bool initialize = false);

    bool getQueue_MS_master(const TlOrbitalInfoObject& orbitalInfo,
                            const bool isCutoffDistribution,
                            const int maxGrainSize,
                            std::vector<Task2>* pTaskList,
                            bool initialize);

    bool getQueue_MS_slave(const TlOrbitalInfoObject& orbitalInfo,
                           const bool isCutoffDistribution,
                           const int maxGrainSize,
                           std::vector<Task2>* pTaskList,
                           bool initialize);

    
    bool getQueue4_DC(const TlOrbitalInfoObject& orbitalInfo,
                      const TlSparseSymmetricMatrix& schwarzTable,
                      const int maxGrainSize,
                      std::vector<Task4>* pTaskList,
                      bool initialize = false);
    bool getQueue4_MS(const TlOrbitalInfoObject& orbitalInfo,
                      const TlSparseSymmetricMatrix& schwarzTable,
                      const int maxGrainSize,
                      std::vector<Task4>* pTaskList,
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

#endif // DFTASKCTRL_PARALLEL_H
