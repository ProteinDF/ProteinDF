#ifndef DFOVERLAP_PARALLEL_H
#define DFOVERLAP_PARALLEL_H

#include "DfOverlap.h"
#include "TlDistributeSymmetricMatrix.h"

class DfOverlap_Parallel : public DfOverlap {
public:
    DfOverlap_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfOverlap_Parallel();

public:
    virtual void getNa(TlVector* pNa);

    virtual void getSpq(TlSymmetricMatrix* pSpq);
    void getSpq(TlDistributeSymmetricMatrix* pSpq);

    virtual void getSab2(TlSymmetricMatrix* pSab2);
    void getSab2(TlDistributeSymmetricMatrix* pSab2);
    
    virtual void getSgd(TlSymmetricMatrix* pSgd);
    void getSgd(TlDistributeSymmetricMatrix* pSgd);
    
    virtual void getdeltaHpqG(const TlVector& deltaMyu, TlSymmetricMatrix& deltaH);
    void getdeltaHpqG(const TlVector& deltaMyu, TlDistributeSymmetricMatrix& deltaH);
    void getdeltaHpqG(const TlDistributeVector& deltaMyu,
                      TlDistributeSymmetricMatrix& deltaH);

protected:
    virtual void logger(const std::string& str) const;
    virtual void cutoffReport();

    virtual std::vector<IJShellPair> getQueue(const TlOrbitalInfo* pOrbitalInfo,
                                              const std::vector<std::vector<std::size_t> >& shellList,
                                              const int maxScore,
                                              const bool initialize = false);

    /// jobの割り当て
    /// 非対称用
    virtual std::vector<IJShellPair> getQueue(const TlOrbitalInfo& orbitalInfo1,
                                              const TlOrbitalInfo& orbitalInfo2,
                                              const ShellListType& shellList1,
                                              const ShellListType& shellList2,
                                              const int maxScore,
                                              const bool initialize = false);
    
    virtual void finalizeIntegral(TlSymmetricMatrix& rMatrix);
    virtual void finalizeIntegral(TlVector& rVector);
};

#endif // DFOVERLAP_PARALLEL_H
