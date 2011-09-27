#ifndef DFJMATRIX_PARALLEL_H
#define DFJMATRIX_PARALLEL_H

#include "DfJMatrix.h"

class DfJMatrix_Parallel : public DfJMatrix {
public:
    DfJMatrix_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfJMatrix_Parallel();

public:
    virtual void buildJMatrix();

protected:
    void buildJMatrix_LAPACK(int iteration,
                             RUN_TYPE runType);
    
    void buildJMatrix_ScaLAPACK(int iteration,
                                RUN_TYPE runType);
};

#endif // DFJMATRIX_PARALLEL_H

