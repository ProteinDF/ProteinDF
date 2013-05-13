#ifndef DFOVERLAPX_PARALLEL_H
#define DFOVERLAPX_PARALLEL_H

#include "DfOverlapX.h"
#include "TlDistributeSymmetricMatrix.h"

class DfOverlapX_Parallel : public DfOverlapX {
public:
    DfOverlapX_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfOverlapX_Parallel();

public:
    void getSpqD(TlDistributeSymmetricMatrix* pSpq);
    void getSabD(TlDistributeSymmetricMatrix* pSab);

    virtual void getSgd(TlSymmetricMatrix* pSgd) {
        DfOverlapX::getSgd(pSgd);
    }
    void getSgd(TlDistributeSymmetricMatrix* pSab);

    /// 変換行列を作成する
    void getTransMat(const TlOrbitalInfoObject& orbInfo1,
                     const TlOrbitalInfoObject& orbInfo2,
                     TlDistributeMatrix* pS);

    /// calc <pq gamma>
    virtual void get_pqg(const TlVector& myu, TlSymmetricMatrix* pF) {
        DfOverlapX::get_pqg(myu, pF);
    }
    virtual void get_pqg(const TlVector& myu, const TlVector& epsilon,
                         TlSymmetricMatrix* pF,
                         TlSymmetricMatrix* pE) {
        DfOverlapX::get_pqg(myu, epsilon, pF, pE);
    }
    virtual void get_pqg(const TlDistributeVector& myu,
                         TlDistributeSymmetricMatrix* pF);
    
    
protected:
    virtual void logger(const std::string& str) const;

    virtual DfTaskCtrl* getDfTaskCtrlObject() const;

    virtual void finalize(TlMatrix* pMtx);
    virtual void finalize(TlSymmetricMatrix* pMtx);
    virtual void finalize(TlVector* pVct);
};

#endif // DFOVERLAPX_PARALLEL_H
