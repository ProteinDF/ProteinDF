#ifndef DFCLEANUP_PARALLEL_H
#define DFCLEANUP_PARALLEL_H

#include "DfCleanup.h"

class DfCleanup_Parallel : public DfCleanup {
public:
    DfCleanup_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfCleanup_Parallel();

    virtual void cleanup();

protected:
    virtual void logger(const std::string& str) const;
};

#endif // DFCLEANUP_PARALLEL_H
