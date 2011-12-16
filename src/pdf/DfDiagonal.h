#ifndef DFDIAGONAL_H
#define DFDIAGONAL_H

#include <string>
#include "DfObject.h"
#include "TlVector.h"

/// KS行列の対角化を行う
class DfDiagonal : public DfObject {
public:
    DfDiagonal(TlSerializeData* pPdfData);
    virtual ~DfDiagonal();

public:
    virtual void DfDiagMain();
    virtual void DfDiagQclo(DfObject::RUN_TYPE runType, const std::string& fragname, int norbcut);

protected:
    template<typename MatrixType, typename SymmetricMatrixType>
    void main(DfObject::RUN_TYPE runType, const std::string& fragname ="", bool bPdfQcloMode = false);

};

// template
template<typename MatrixType, typename SymmetricMatrixType>
void DfDiagonal::main(const DfObject::RUN_TYPE runType, const std::string& fragname, bool bPdfQcloMode)
{
    // input file
    SymmetricMatrixType Fprime = this->getFprimeMatrix<SymmetricMatrixType>(runType, this->m_nIteration, fragname);

    // diagonal
    TlVector* pEigVal = new TlVector();
    MatrixType* pEigVec = new MatrixType();
    Fprime.diagonal(pEigVal, pEigVec);

    // save eigval
    pEigVal->save(this->getEigenvaluesPath(runType, this->m_nIteration));

    // save eigvec
    this->saveCprimeMatrix(runType, this->m_nIteration, fragname, *pEigVec);

    // finalize
    delete pEigVal;
    pEigVal = NULL;

    delete pEigVec;
    pEigVec = NULL;
}

#endif // DFDIAGONAL_H

