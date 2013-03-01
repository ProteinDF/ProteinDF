#ifndef DFINVMATRIX_PARALLEL_H
#define DFINVMATRIX_PARALLEL_H

#include "DfInvMatrix.h"

class DfInvMatrix_Parallel : public DfInvMatrix {
public:
    DfInvMatrix_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfInvMatrix_Parallel();

    virtual void DfInvMain();
};

#endif // DFINVMATRIX_PARALLEL_H
