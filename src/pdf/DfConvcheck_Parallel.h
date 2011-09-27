#ifndef DFCONVCHECK_PARALLEL_H
#define DFCONVCHECK_PARALLEL_H

#include "DfConvcheck.h"

class DfConvcheck_Parallel : public DfConvcheck {
public:
    DfConvcheck_Parallel(TlSerializeData* pPdfParam, int num_iter);
    virtual ~DfConvcheck_Parallel();

public:
    virtual void DfConvcheckMain();

protected:
    virtual void logger(const std::string& str) const;
    void main_ScaLAPACK();
};

#endif // DFCONVCHECK_PARALLEL_H
