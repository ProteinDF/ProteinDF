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
    void getK_conventional(const RUN_TYPE runType, TlSymmetricMatrix *pK);
    void getK_byCD(const RUN_TYPE runType, TlSymmetricMatrix *pK);
    void getK_byRI_K(const RUN_TYPE runType, TlSymmetricMatrix *pK);
};

#endif // DFKMATRIX_H
