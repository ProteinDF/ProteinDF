#ifndef DFJMATRIX_H
#define DFJMATRIX_H

#include "DfObject.h"
#include "TlSymmetricMatrix.h"
#include "TlVector.h"

class DfJMatrix : public DfObject {
public:
    DfJMatrix(TlSerializeData* pPdfParam);
    virtual ~DfJMatrix();

public:
    virtual void buildJ();

protected:
    virtual void getJ_RI();
    virtual void getJ_CD();
    virtual void getJ_conventional();

    virtual void getJ_RI_local(TlSymmetricMatrix* pJ);
    virtual void getJ_CD_local(TlSymmetricMatrix* pJ);
    virtual void getJ_conventional_local(TlSymmetricMatrix* pJ);

protected:
    virtual void saveJMatrix(const TlSymmetricMatrix& J);
    virtual TlSymmetricMatrix getJMatrix(const int iteration);

    virtual TlVector getRho(const RUN_TYPE runType, const int iteration);
    virtual TlSymmetricMatrix getPMatrix(const int iteration);
    virtual TlSymmetricMatrix getDiffDensityMatrix();
};

#endif // DFJMATRIX_H
