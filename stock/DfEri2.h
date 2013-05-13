#ifndef DFERI2_H
#define DFERI2_H

#include "DfObject.h"
#include "TlSymmetricMatrix.h"

class DfEri;

class DfEri2 : public DfObject {
public:
    DfEri2(TlSerializeData* pPdfParam);
    virtual ~DfEri2();

public:
    TlSymmetricMatrix generateInvSquareVMatrix();
    TlSymmetricMatrix getKMatrix();
    TlSymmetricMatrix getKMatrix(const TlSymmetricMatrix& P);

protected:
    virtual DfEri* getDfEri();
    virtual TlSymmetricMatrix loadInvSquareVMatrix();
    virtual TlMatrix loadCutoffCMatrix();

protected:
    // to implement for parallel
    virtual TlSymmetricMatrix getPMatrix(int iteration);

    TlSymmetricMatrix getKMatrixUsingCMatrix(const TlMatrix& C);
    TlSymmetricMatrix getKMatrixUsingDensityMatrix(const TlSymmetricMatrix& P);

    virtual TlSymmetricMatrix loadKMatrix(const int nIteration);
    void saveKMatrix(const TlSymmetricMatrix& K, const int nIteration);

private:
    TlSymmetricMatrix m_transformedIntegral;

    bool m_bUseDensityMatrix;
};

#endif // DFERI2_H

