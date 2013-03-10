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
    virtual TlSymmetricMatrix getKMatrix(const RUN_TYPE runType,
                                         const int iteration);
    virtual void saveKMatrix(const RUN_TYPE runType,
                             const TlSymmetricMatrix& K);

    template<class SymmetricMatrixType>
    SymmetricMatrixType getDiffDensityMatrix(RUN_TYPE runType);

    template<class SymmetricMatrixType>
    SymmetricMatrixType getDensityMatrix(RUN_TYPE runType);
};

template<class SymmetricMatrixType>
SymmetricMatrixType DfKMatrix::getDiffDensityMatrix(const RUN_TYPE runType)
{
    SymmetricMatrixType diffP;
    switch (runType) {
    case RUN_RKS:
        diffP = DfObject::getDiffDensityMatrix<SymmetricMatrixType>(RUN_RKS, this->m_nIteration);
        diffP *= 0.5;
        break;
        
    default:
        // RUN_UKS_ALPHA, RUN_UKS_BETA, RUN_ROKS_CLOSE, RUN_ROKS_OPEN
        diffP  = DfObject::getDiffDensityMatrix<SymmetricMatrixType>(runType, this->m_nIteration);
        break;
    }
        
    return diffP;
}

template<class SymmetricMatrixType>
SymmetricMatrixType DfKMatrix::getDensityMatrix(const RUN_TYPE runType)
{
    SymmetricMatrixType P;
    switch (runType) {
    case RUN_RKS:
        P = DfObject::getPpqMatrix<SymmetricMatrixType>(RUN_RKS, this->m_nIteration -1);
        P *= 0.5;
        break;
        
    default:
        // RUN_UKS_ALPHA, RUN_UKS_BETA, RUN_ROKS_CLOSE, RUN_ROKS_OPEN
        P = DfObject::getPpqMatrix<SymmetricMatrixType>(runType, this->m_nIteration -1);
        break;
    }
        
    return P;
}

#endif // DFKMATRIX_H
