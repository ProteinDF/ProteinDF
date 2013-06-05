#ifndef DFERI_PARALLEL_H
#define DFERI_PARALLEL_H

#include "DfEri.h"
#include "TlDistributeVector.h"
#include "TlDistributeSymmetricMatrix.h"

class DfEri_Parallel : public DfEri {
public:
    DfEri_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfEri_Parallel();

public:
    virtual void getSab(TlSymmetricMatrix* pSab);
    void getSab(TlDistributeSymmetricMatrix* pSab);

    virtual void getDeltaT(const TlSymmetricMatrix& rDeltaPpq, TlVector* pDeltaT);
    void getDeltaT(const TlDistributeSymmetricMatrix& rDeltaPpq, TlDistributeVector* pDeltaT);

    virtual void getdeltaHpqA(const TlVector& deltaRho, TlSymmetricMatrix& deltaH);
    void getdeltaHpqA(const TlVector& deltaRho,
                      TlDistributeSymmetricMatrix& deltaH);
    void getdeltaHpqA(const TlDistributeVector& deltaRho,
                      TlDistributeSymmetricMatrix& deltaH);

protected:
    void getDeltaT_core(const TlDistributeSymmetricMatrix& deltaPpq,
                        TlDistributeVector* pDeltaT);

    void getDeltaHpqA_rev2(const TlVector& deltaRho,
                           TlDistributeSymmetricMatrix& deltaH);

    
    virtual std::vector<IJShellPair> getQueue(const TlOrbitalInfo* pOrbitalInfo,
                                              const std::vector<std::vector<std::size_t> >& shellList,
                                              const int maxNumOfPairs, const bool initialize = false);

    virtual void finalizeIntegral(TlSymmetricMatrix& rMatrix);
    virtual void finalizeIntegral(TlVector& rVector);

    virtual void finalizeIntegral(TlDistributeSymmetricMatrix& rMatrix) {
        // do nothing
    }

protected:
    virtual void cutoffReport();

protected:
    // tag for MPI
    enum {
        TAG_REQUEST_JOB = 991,
        TAG_ASSIGN_JOB = 992,
        TAG_TERMINATE_SLAVE = 993,
        TAG_TERMINATE_OK = 994
    };
    
    struct DT_Job {
        DT_Job(const TlSparseSymmetricMatrix& p,
               const std::vector<IJShellPair>& ijs)
            : partP(p), ijShellPairs(ijs) {
        }
        TlSparseSymmetricMatrix partP;
        std::vector<IJShellPair> ijShellPairs;
    };
    
protected:
    std::size_t densityMatrixCacheMemSize_;

    /// 1プロセスあたりのタスク数
    ///
    /// 分散行列を用いたmaster-slave型計算において、同時に通信待ちするタスクの数
    int numOfTasksPerProc_;

    bool isDebugOut_getDeltaT_;
};

#endif // DFERI_PARALLEL_H
