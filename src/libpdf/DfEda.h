#ifndef DFEDA_H
#define DFEDA_H

#include "DfObject.h"

class DfGridFreeXC;

class DfEda : public DfObject {
public:
    DfEda(TlSerializeData* pPdfParam);
    virtual ~DfEda();

public:
    void calc();

protected:
    void calcGridFreeXC();
    virtual DfGridFreeXC* getDfGridFreeXcObject();
};

#endif  // DFEDA_H
