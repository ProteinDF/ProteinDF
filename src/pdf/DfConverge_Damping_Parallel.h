#ifndef DFCONVERGE_DAMPING_PARALLEL_H
#define DFCONVERGE_DAMPING_PARALLEL_H

#include "DfConverge_Damping.h"

class DfConverge_Damping_Parallel : public DfConverge_Damping {
public:
    DfConverge_Damping_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfConverge_Damping_Parallel();

protected:
    virtual void logger(const std::string& str) const;
    
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


#endif // DFCONVERGE_DAMPING_PARALLEL_H
