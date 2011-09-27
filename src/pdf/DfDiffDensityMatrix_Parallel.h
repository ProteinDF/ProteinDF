#ifndef DFDIFFDENSITYMATRIX_PARALLEL_H
#define DFDIFFDENSITYMATRIX_PARALLEL_H

#include "DfDiffDensityMatrix.h"

class DfDiffDensityMatrix_Parallel : public DfDiffDensityMatrix {
public:
    DfDiffDensityMatrix_Parallel(TlSerializeData* pPdfParam);
    ~DfDiffDensityMatrix_Parallel();

public:
    virtual void exec();

protected:
    virtual void logger(const std::string& str) const;
    void calc_usingScaLAPACK(DfObject::RUN_TYPE runType, int iteration);
};

#endif // DFDIFFDENSITYMATRIX_PARALLEL_H
