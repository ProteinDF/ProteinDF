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

    template<class SymmetricMatrixType>
    SymmetricMatrixType getDiffDensityMatrix();

    template<class SymmetricMatrixType>
    SymmetricMatrixType getDensityMatrix();
};

template<class SymmetricMatrixType>
SymmetricMatrixType DfJMatrix::getDiffDensityMatrix()
{
    SymmetricMatrixType diffP;
    switch (this->m_nMethodType) {
    case METHOD_RKS:
        diffP = DfObject::getDiffDensityMatrix<SymmetricMatrixType>(RUN_RKS, this->m_nIteration);
        break;
        
    case METHOD_UKS:
        diffP  = DfObject::getDiffDensityMatrix<SymmetricMatrixType>(RUN_UKS_ALPHA, this->m_nIteration);
        diffP += DfObject::getDiffDensityMatrix<SymmetricMatrixType>(RUN_UKS_BETA,  this->m_nIteration);
        break;

    case METHOD_ROKS:
        diffP  = DfObject::getDiffDensityMatrix<SymmetricMatrixType>(RUN_ROKS_CLOSE, this->m_nIteration);
        diffP += DfObject::getDiffDensityMatrix<SymmetricMatrixType>(RUN_ROKS_OPEN,  this->m_nIteration);
        break;

    default:
        this->log_.critical("program error");
        break;
    }
        
    return diffP;
}

template<class SymmetricMatrixType>
SymmetricMatrixType DfJMatrix::getDensityMatrix()
{
    SymmetricMatrixType P;
    switch (this->m_nMethodType) {
    case METHOD_RKS:
        P = DfObject::getPpqMatrix<SymmetricMatrixType>(RUN_RKS, this->m_nIteration -1);
        break;
        
    case METHOD_UKS:
        P  = DfObject::getPpqMatrix<SymmetricMatrixType>(RUN_UKS_ALPHA, this->m_nIteration -1);
        P += DfObject::getPpqMatrix<SymmetricMatrixType>(RUN_UKS_BETA,  this->m_nIteration -1);
        break;

    case METHOD_ROKS:
        P  = DfObject::getPpqMatrix<SymmetricMatrixType>(RUN_ROKS_CLOSE, this->m_nIteration -1);
        P += DfObject::getPpqMatrix<SymmetricMatrixType>(RUN_ROKS_OPEN,  this->m_nIteration -1);
        break;

    default:
        this->log_.critical("program error");
        break;
    }
        
    return P;
}

#endif // DFJMATRIX_H
