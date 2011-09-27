#ifndef DFHPQX_PARALLEL_H
#define DFHPQX_PARALLEL_H

#include "DfHpqX.h"
#include "TlDistributeSymmetricMatrix.h"

class DfHpqX_Parallel : public DfHpqX {
public:
    DfHpqX_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfHpqX_Parallel();

public:
    void getHpqD(TlDistributeSymmetricMatrix* pHpq,
                 TlDistributeSymmetricMatrix* pHpq2);
    
protected:
    virtual void logger(const std::string& str) const;

    virtual DfTaskCtrl* getDfTaskCtrlObject() const;

    virtual void finalize(TlSymmetricMatrix* pHpq, TlSymmetricMatrix* pHpq2);
};

#endif // DFHPQX_PARALLEL_H
