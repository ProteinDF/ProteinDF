#ifndef DFJMATRIX_H
#define DFJMATRIX_H

#include "DfObject.h"

class DfJMatrix : public DfObject {
public:
    DfJMatrix(TlSerializeData* pPdfParam);
    virtual ~DfJMatrix();

public:
    virtual void buildJMatrix();
};

#endif // DFJMATRIX_H
