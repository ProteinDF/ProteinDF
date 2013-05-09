#ifndef DFXMATRIX_PARALLEL_H
#define DFXMATRIX_PARALLEL_H

#include "DfXMatrix.h"

class DfXMatrix_Parallel : public DfXMatrix {
public:
    DfXMatrix_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfXMatrix_Parallel();

public:
    virtual void buildX();

protected:
    void buildX_LAPACK();
    void buildX_ScaLAPACK();
};

#endif // DFXMATRIX_PARALLEL_H
