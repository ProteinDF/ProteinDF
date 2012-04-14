#ifndef DFJMATRIX_PARALLEL_H
#define DFJMATRIX_PARALLEL_H

#include "DfJMatrix.h"

class DfJMatrix_Parallel : public DfJMatrix {
public:
    DfJMatrix_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfJMatrix_Parallel();

protected:
    virtual void getJ_RI();
    virtual void getJ_CD();
    virtual void getJ_conventional();

    virtual void getJ_RI_local(TlSymmetricMatrix* pJ);
    void getJ_RI_distributed(TlDistributeSymmetricMatrix *pJ);

    virtual void getJ_CD_local(TlSymmetricMatrix* pJ);
    void getJ_CD_distributed(TlDistributeSymmetricMatrix *pJ);

    virtual void getJ_conventional_local(TlSymmetricMatrix *pJ);
    void getJ_conventional_distributed(TlDistributeSymmetricMatrix *pJ);

protected:
    virtual void saveJMatrix(const TlSymmetricMatrix& J);
    virtual TlSymmetricMatrix getJMatrix(const int iteration);

    virtual TlVector getRho(const RUN_TYPE runType, const int iteration);
    virtual TlSymmetricMatrix getPMatrix(const int iteration);
    virtual TlSymmetricMatrix getDiffDensityMatrix();

};

#endif // DFJMATRIX_PARALLEL_H

