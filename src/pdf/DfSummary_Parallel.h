#ifndef DFSUMMARY_PARALLEL_H
#define DFSUMMARY_PARALLEL_H

#include "DfSummary.h"

class DfSummary_Parallel : public DfSummary {
public:
    DfSummary_Parallel(TlSerializeData* pPdfParam_);
    ~DfSummary_Parallel();

    void exec();

protected:
    virtual void logger(const std::string& str) const;
    
    void exec_LAPACK();
    void exec_ScaLAPACK();
};

#endif // DFSUMMARY_PARALLEL_H
