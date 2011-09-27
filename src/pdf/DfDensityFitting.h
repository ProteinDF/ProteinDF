#ifndef DFDENSITYFITTING_H
#define DFDENSITYFITTING_H

#include "DfObject.h"
#include "TlVector.h"
#include "TlSymmetricMatrix.h"
#include "DfDensityFittingObject.h"
#include "DfEri.h"
#include "DfEriX.h"

class DfDensityFitting : public DfDensityFittingTmpl<TlSymmetricMatrix, TlVector, DfEri> {
public:
    DfDensityFitting(TlSerializeData* pPdfParam);
    virtual ~DfDensityFitting();

public:
    void exec();
};


class DfDensityFittingX : public DfDensityFittingTmpl<TlSymmetricMatrix, TlVector, DfEriX> {
public:
    DfDensityFittingX(TlSerializeData* pPdfParam);
    virtual ~DfDensityFittingX();

public:
    void exec();
};

#endif // DFDENSITYFITTING_H
