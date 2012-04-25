#ifndef DFCLEANUP_PARALLEL_H
#define DFCLEANUP_PARALLEL_H

#include "DfCleanup.h"

class DfCleanup_Parallel : public DfCleanup {
public:
    DfCleanup_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfCleanup_Parallel();

    virtual void cleanup();
};

#endif // DFCLEANUP_PARALLEL_H
