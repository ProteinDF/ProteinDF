#ifndef DFEDA_PARALLEL_H
#define DFEDA_PARALLEL_H

#include "DfEda.h"

class DfEda_Parallel : public DfEda {
public:
    DfEda_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfEda_Parallel();

protected:
    virtual DfGridFreeXC* getDfGridFreeXcObject();
};

#endif  // DFEDA_PARALLEL_H
