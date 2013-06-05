#ifndef DFHPQ_PARALLEL_H
#define DFHPQ_PARALLEL_H

#include "DfHpq.h"
#include "TlSymmetricMatrix.h"

class DfHpq_Parallel : public DfHpq {
public:
    DfHpq_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfHpq_Parallel();

public:
    /// Hpqを算出する
    /// 
    /// @param [out] pHoq Hpq行列を受け取るポインタ
    /// @param [out] pHoq2 Hpq2行列を受け取るポインタ
    void getHpq(TlSymmetricMatrix* pHpq,
                       TlSymmetricMatrix* pHpq2);

    /// 大域分散行列でHpq, Hp2を算出する
    ///
    /// @param [out] pHoq Hpq行列を受け取るポインタ
    /// @param [out] pHoq2 Hpq2行列を受け取るポインタ
    void getHpq(TlDistributeSymmetricMatrix* pHpq,
                TlDistributeSymmetricMatrix* pHpq2);
    
protected:
    virtual std::vector<IJShellPair> getQueue(int maxNumOfPairs, bool initialize = false);
    std::vector<IJShellPair> getQueue_DC(int maxNumOfPairs, bool initialize);
    std::vector<IJShellPair> getQueue_MS(int maxNumOfPairs, bool initialize);

protected:
    virtual void logger(const std::string& str) const;
    virtual void parallelLogger(const std::string& str);

protected:
    enum {
        TAG_HPQ_REQUEST_TASK = 1001,
        TAG_HPQ_SEND_TASK = 1002
    };
    
protected:
    int MS_blockSize_;
};

#endif // DFHPQ_PARALLEL_H
