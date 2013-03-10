#ifndef DFKMATRIX_PARALLEL_H
#define DFKMATRIX_PARALLEL_H

#include "DfKMatrix.h"

class DfKMatrix_Parallel : public DfKMatrix {
public:
    DfKMatrix_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfKMatrix_Parallel();

protected:
    // virtual void getK_RI();
    virtual void getK_CD();
    virtual void getK_conventional();

    // virtual void getK_RI_local(const RUN_TYPE runType, TlSymmetricMatrix* pJ);
    // void getK_RI_distributed(const RUN_TYPE runType, TlDistributeSymmetricMatrix *pJ);

    virtual void getK_CD_local(const RUN_TYPE runType, TlSymmetricMatrix* pJ);
    void getK_CD_distributed(const RUN_TYPE runType, TlDistributeSymmetricMatrix *pJ);

    virtual void getK_conventional_local(const RUN_TYPE runType, TlSymmetricMatrix *pJ);
    void getK_conventional_distributed(const RUN_TYPE runType, TlDistributeSymmetricMatrix *pJ);

protected:
    // virtual TlSymmetricMatrix getPMatrix(const RUN_TYPE runType,
    //                                      const int iteration);
    // virtual TlSymmetricMatrix getDiffDensityMatrix(const RUN_TYPE runType);

    virtual TlSymmetricMatrix getKMatrix(const RUN_TYPE runType,
                                         const int iteration);
    virtual void saveKMatrix(const RUN_TYPE runType,
                             const TlSymmetricMatrix& K);
};

#endif // DFKMATRIX_PARALLEL_H

