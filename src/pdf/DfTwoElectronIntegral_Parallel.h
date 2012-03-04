#ifndef DFTWOELECTRONINTEGRAL_PARALLEL_H
#define DFTWOELECTRONINTEGRAL_PARALLEL_H

#ifdef USE_OLD_ERI_ENGINE

#include "DfTwoElectronIntegral.h"
#include "TlDistributeSymmetricMatrix.h"
#include "TlSparseVectorMatrix.h"

class DfTwoElectronIntegral_Parallel : public DfTwoElectronIntegral {
public:
    DfTwoElectronIntegral_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfTwoElectronIntegral_Parallel();

    virtual void getContractKMatrixByIntegralDriven(const TlSymmetricMatrix& P,
                                                    TlSymmetricMatrix* pK);
    virtual void getContractKMatrixByRTmethod(const TlSymmetricMatrix& P,
                                              TlSymmetricMatrix* pK);

    void getContractKMatrixByRTmethod(const TlDistributeSymmetricMatrix& P,
                                      TlDistributeSymmetricMatrix* pK);
    
protected:
    virtual void logger(const std::string& str) const;
    void loggerP(const std::string& str) const;
    
    // this->m_IKShellPairListから計算すべきIK shellの組を返す
    // 並列化用(シリアルでは全部返す)
    virtual std::vector<DfTwoElectronIntegral::IKShellPair> getLocalIKShellPairList(const int nIKShellPairIndex);

    virtual void finalize(TlSymmetricMatrix& K);

protected:
    // for Divide & Conquire


protected:
    // for Master-Slave
    void getContractKMatrixByIntegralDrivenCore_MS(const TlSymmetricMatrix& P,
                                                   TlSymmetricMatrix* pK);
    bool assignJob(unsigned int* pEriType,
                   int* pIShell, int* pKShell,
                   double* pProgress = NULL,
                   bool isInitialize = false);

protected:
    // ScaLAPACK ===============================================================
    // メッセージコマンド
    enum {
        REQUEST_DENS_VECTOR = -1,
        REQUEST_JOB = -2,
        TERMINATE_JOB = -3,
        COMPLETION_JOB = -4
    };

    // MPI tags
    enum {
        TAG_COMMAND = 1,
        TAG_REQUEST_VECTOR = 2,
        TAG_VECTOR_CONTENTS = 3
    };
    
    
    // 計算タスク用構造体
    struct CalcShellPair {
    public:
        CalcShellPair(int sID, const TlOrbitalInfo& orbInfo,
                      const index_type numOfAOs,
                      unsigned int et,
                      index_type i, index_type k)
            : sessionId(sID),
              Pi(numOfAOs), Pk(numOfAOs),
              eriType(et), i_shell(i), k_shell(k),
              isReadyForPi(false), isReadyForPk(false) {

            const index_type startI = this->i_shell;
            const index_type endI = startI + (orbInfo.getShellType(i) * 2 + 1);
            for (index_type I = startI; I < endI; ++I) {
                for (index_type col = 0; col < numOfAOs; ++col) {
                    this->Pi.set(I, col, 0.0);
                }
            }
            
            const index_type startK = this->k_shell;
            const index_type endK = startK + (orbInfo.getShellType(k) * 2 + 1);
            for (index_type K = startK; K < endK; ++K) {
                for (index_type col = 0; col < numOfAOs; ++col) {
                    this->Pk.set(K, col, 0.0);
                }
            }
        }

    public:
        int sessionId;
        TlSparseSymmetricMatrix Pi;
        TlSparseSymmetricMatrix Pk;
        unsigned int eriType;
        int i_shell;
        int k_shell;
        bool isReadyForPi;
        bool isReadyForPk;
    };
    typedef std::list<CalcShellPair> CalcTasksType;

    // ベクトル送信用構造体
    struct SendVectorInfo {
    public:
        SendVectorInfo(int t, int s) : to(t), shell(s), pBuf(NULL) {
        }
        
    public:
        int to;
        int shell;
        double* pBuf;
    };
    typedef std::list<SendVectorInfo> SendVectorTasksType;

    void RT_master(TlSparseVectorMatrix* pP,
                   const std::map<int, int>& densityVectorHoldingTable);
    void RT_slave(TlSparseVectorMatrix* pP,
                  TlSparseSymmetricMatrix* pK,
                  const std::map<int, int>& densityVectorHoldingTable);

    void RT_master2(const TlDistributeSymmetricMatrix& P);
    void RT_slave2(const TlDistributeSymmetricMatrix& P,
                   TlSparseSymmetricMatrix* pK);


    
    void requestDensityVector(const std::map<int, int>& densityVectorHoldingTable,
                              const int shellBegin,
                              std::map<int, int*>* pRequestDensVecCmdTbl,
                              std::map<int, double*>* pRecvDensVecCmdTbl);
    void cleanupRequestDensVecCmdTbl(std::map<int, int*>* pRequestDensVecCmdTbl);
    void cleanupRecvDensVecCmdTbl(std::map<int, double*>* pRecvDensVecCmdTbl,
                                  TlSparseVectorMatrix* pP);

    void sendDensityVector(const TlSparseVectorMatrix& P,
                           SendVectorTasksType* pSendVectorTasks);

protected:
    enum {
        TAG_TEI_REQUEST_TASK = 1101,
        TAG_TEI_ASSIGN_TASK = 1102
    };
    
protected:
    std::size_t densityMatrixCacheMemSize_;

    /// 2電子積分オリジナルの通信方法を選択
    bool isTeiOriginalCommunication_;

    /// 1プロセスあたりの保持するタスクの量
    int numOfTasksPerProc_;

    /// デバッグ出力のスイッチ
    bool isDebugOutMsg_;
};

#endif // USE_OLD_ERI_ENGINE
#endif // DFTWOELECTRONINTEGRAL_PARALLEL_H
