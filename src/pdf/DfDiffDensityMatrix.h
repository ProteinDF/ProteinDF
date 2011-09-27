#ifndef DFDIFFDENSITYMATRIX_H
#define DFDIFFDENSITYMATRIX_H

#include "DfObject.h"
#include "TlSymmetricMatrix.h"

class DfDiffDensityMatrix : public DfObject {
public:
    DfDiffDensityMatrix(TlSerializeData* pPdfParam);
    virtual ~DfDiffDensityMatrix();

public:
    /// 差電子密度行列を求める
    ///
    virtual void exec();

protected:
    void calc(DfObject::RUN_TYPE runType, int iteration);

    template<class SymmetricMatrixType>
    void calc_ROKS();

protected:
    /// 差電子密度行列をディスクに保存する(true)かどうか
    bool isSaveDiffMatrix_;
};


template<class SymmetricMatrixType>
void DfDiffDensityMatrix::calc_ROKS()
{
    SymmetricMatrixType P_close;
    P_close.load(DfObject::getP1pqMatrixPath(this->m_nIteration -1));
    SymmetricMatrixType P_open;
    P_open.load(DfObject::getP2pqMatrixPath(this->m_nIteration -1));
    P_open += P_close;

    {
        SymmetricMatrixType prevP_open;
        prevP_open.load(DfObject::getP2pqMatrixPath(this->m_nIteration -2));
        P_open -= prevP_open;
    }

    {
        SymmetricMatrixType prevP_close;
        prevP_close.load(DfObject::getP1pqMatrixPath(this->m_nIteration -2));
        P_open -= prevP_close;
        P_close -= prevP_close;
    }

    P_open.save(this->getDiffDensityMatrixPath(DfObject::RUN_ROKS_OPEN, this->m_nIteration));
    P_close.save(this->getDiffDensityMatrixPath(DfObject::RUN_ROKS_CLOSE, this->m_nIteration));
}


#endif // DFDIFFDENSITYMATRIX_H
