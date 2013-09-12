#ifndef DFXMATRIX_PARALLEL_H
#define DFXMATRIX_PARALLEL_H

#include "DfXMatrix.h"

class DfXMatrix_Parallel : public DfXMatrix {
public:
    DfXMatrix_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfXMatrix_Parallel();

public:
    virtual void buildX();
    
    virtual void canonicalOrthogonalize(const TlSymmetricMatrix& S,
                                        TlMatrix* pX, TlMatrix* pXinv);

    virtual void lowdinOrthogonalize(const TlSymmetricMatrix& S,
                                     TlMatrix* pX, TlMatrix* pXinv);

    void canonicalOrthogonalize(const TlDistributeSymmetricMatrix& S,
                                TlDistributeMatrix* pX,
                                TlDistributeMatrix* pXinv);

    void lowdinOrthogonalize(const TlDistributeSymmetricMatrix& S,
                             TlDistributeMatrix* pX,
                             TlDistributeMatrix* pXinv);

protected:
    void buildX_LAPACK();
    void buildX_ScaLAPACK();
};

#endif // DFXMATRIX_PARALLEL_H
