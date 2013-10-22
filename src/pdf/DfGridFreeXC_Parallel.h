#ifndef DFGRIDFREEXC_PARALLEL_H
#define DFGRIDFREEXC_PARALLEL_H

#include "DfGridFreeXC.h"

class DfGridFreeXC_Parallel : public DfGridFreeXC 
{
public:
    DfGridFreeXC_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfGridFreeXC_Parallel();

public:
    virtual void preprocessBeforeSCF();
    
protected:
    virtual DfOverlapX* getDfOverlapObject();
    virtual DfXMatrix* getDfXMatrixObject();

    void preprocessBeforeSCF_LAPACK();
    void preprocessBeforeSCF_ScaLAPACK();

protected:
    virtual void buildFxc_LDA();
    void buildFxc_LDA_LAPACK();
    void buildFxc_LDA_ScaLAPACK();
};

#endif // DFGRIDFREEXC_PARALLEL_H
