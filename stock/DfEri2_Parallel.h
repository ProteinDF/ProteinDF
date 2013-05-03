#ifndef DFERI2_PARALLEL_H
#define DFERI2_PARALLEL_H

#include "DfEri2.h"

class DfEri2_Parallel : public DfEri2 {
public:
    DfEri2_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfEri2_Parallel();

protected:
    virtual DfEri* getDfEri();
    virtual TlSymmetricMatrix loadInvSquareVMatrix();
    virtual TlMatrix loadCutoffCMatrix();

protected:
    // to implement for parallel
    virtual TlSymmetricMatrix getPMatrix(int iteration);

    //TlSymmetricMatrix getKMatrixUsingCMatrix(const TlMatrix& C);
    //TlSymmetricMatrix getKMatrixUsingDensityMatrix(const TlSymmetricMatrix& P);

    virtual TlSymmetricMatrix loadKMatrix(const int nIteration);
    virtual void saveKMatrix(const TlSymmetricMatrix& K, const int nIteration);

};

#endif // DFERI2_PARALLEL_H

