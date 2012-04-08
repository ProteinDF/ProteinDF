#ifndef DFJMATRIX_H
#define DFJMATRIX_H

#include "DfObject.h"

class DfJMatrix : public DfObject {
public:
    DfJMatrix(TlSerializeData* pPdfParam);
    virtual ~DfJMatrix();

public:
    virtual void buildJ();

protected:
    void getJ_conventional(TlSymmetricMatrix* pJ);
    void getJ_RI(TlSymmetricMatrix* pJ);
    void getJ_CD(TlSymmetricMatrix* pJ);
};

#endif // DFJMATRIX_H
