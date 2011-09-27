#ifndef DFLOCALIZE_PARALLEL_H
#define DFLOCALIZE_PARALLEL_H

#include <vector>
#include "DfLocalize.h"

class DfLocalize_Parallel : public DfLocalize {
public:
    DfLocalize_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfLocalize_Parallel();

    virtual void localize();

protected:
    virtual int getJobItem(DfLocalize::JobItem* pJob, bool isInitialized = false);

protected:
    std::vector<bool> jobFinishedList_;
    std::vector<bool> jobOccupiedOrb_;
};

#endif // DFLOCALIZE_PARALLEL_H
