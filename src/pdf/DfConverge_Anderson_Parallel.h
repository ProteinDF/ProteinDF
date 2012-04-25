#ifndef DFCONVERGE_ANDERSON_PARALLEL_H
#define DFCONVERGE_ANDERSON_PARALLEL_H

#include "DfConverge_Anderson.h"

class DfConverge_Anderson_Parallel : public DfConverge_Anderson {
public:
    DfConverge_Anderson_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfConverge_Anderson_Parallel();

protected:
    virtual void convergeRhoTilde();
    virtual void convergeKSMatrix();
    virtual void convergePMatrix();

    void convergeRhoTilde_LAPACK();
    void convergeKSMatrix_LAPACK();
    void convergePMatrix_LAPACK();

    void convergeRhoTilde_ScaLAPACK();
    void convergeKSMatrix_ScaLAPACK();
    void convergePMatrix_ScaLAPACK();
};

#endif // DFCONVERGE_ANDERSON_PARALLEL_H
