#ifndef DFERIOBJECT_H
#define DFERIOBJECT_H

#include "DfObject.h"

class DfEriObject : public DfObject {
    DfEriObject(TlSerializeData* pPdfParam);
    virtual ~DfEriObject();

public:
    virtual void getJ(const TlSymmetricMatrix& P, TlVector* pRho);
    

};

#endif // DFERIOBJECT_H
