#ifndef DFCONVERGE_H
#define DFCONVERGE_H

#include "DfObject.h"
#include "TlVector.h"
#include "TlSymmetricMatrix.h"

class DfConverge : public DfObject {
protected:
    enum ConvergeTarget {
        RHO_TILDE,
        KS_MATRIX,
        DENSITY_MATRIX
    };

public:
    DfConverge(TlSerializeData* pPdfParam);
    virtual ~DfConverge();

public:
    void doConverge();

protected:
    virtual void convergeRhoTilde() =0;
    virtual void convergeKSMatrix() =0;
    virtual void convergePMatrix() =0;

protected:
    ConvergeTarget m_nConvergeTarget;
};

#endif // DFCONVERGE
