#ifndef DFKMATRIX_H
#define DFKMATRIX_H

#include "DfObject.h"

class DfKMatrix : public DfObject {
public:
    DfKMatrix(TlSerializeData* pPdfParam);
    virtual ~DfKMatrix();

public:
    void buildK();

protected:
    //virtual void getK_RI();
    virtual void getK_CD();
    virtual void getK_conventional();

    //void getK_RI_local(const RUN_TYPE runType, TlSymmetricMatrix *pK);
    void getK_CD_local(const RUN_TYPE runType, TlSymmetricMatrix *pK);
    void getK_conventional_local(const RUN_TYPE runType, TlSymmetricMatrix *pK);

protected:
    virtual TlSymmetricMatrix getPMatrix(const RUN_TYPE runType,
                                         const int iteration);
    virtual TlSymmetricMatrix getDiffDensityMatrix(const RUN_TYPE runType);
    virtual TlSymmetricMatrix getKMatrix(const RUN_TYPE runType,
                                         const int iteration);
    virtual void saveKMatrix(const RUN_TYPE runType,
                             const TlSymmetricMatrix& K);

};

#endif // DFKMATRIX_H
