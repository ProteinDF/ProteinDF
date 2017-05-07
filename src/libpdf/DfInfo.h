#ifndef DFINFO_H
#define DFINFO_H

#include "DfObject.h"

class DfInfo : public DfObject {
public:
    DfInfo(TlSerializeData* pPdfParam);
    virtual ~DfInfo();

public:
    int getNumOfElectrons() const;
};

#endif // DFINFO_H
