#ifndef DFDENSITYFITTINGX_H
#define DFDENSITYFITTINGX_H

#include "DfObject.h"
#include "TlVector.h"
#include "TlSymmetricMatrix.h"
#include "DfDensityFittingObject.h"
#include "DfEriX.h"

class DfDensityFittingX : public DfDensityFittingTmpl<TlSymmetricMatrix, TlVector, DfEriX> {
public:
    DfDensityFittingX(TlSerializeData* pPdfParam);
    virtual ~DfDensityFittingX();

public:
    void exec();
};

#endif // DFDENSITYFITTINGX_H
