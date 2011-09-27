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
    
protected:
    virtual void logger(const std::string& str) const;

    virtual DfTaskCtrl* getDfTaskCtrlObject() const;

    virtual void finalize(TlMatrix* pMtx);
    virtual void finalize(TlSymmetricMatrix* pMtx);
    virtual void finalize(TlVector* pVct);
};

#endif // DFOVERLAPX_PARALLEL_H
