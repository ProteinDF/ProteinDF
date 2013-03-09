#ifndef DFINVMATRIX_H
#define DFINVMATRIX_H

#include <string>
#include "DfObject.h"

// Function : Calculate inverse matrix by means of modified Cholesky decomposion (MCD)
//            first, evaluate inverse matrix of cd basis overlap matrix (fl_Mtr_Sabinv)
//            next, evaluate inverse matrix of xc fitting basis overlap matrix (fl_Mtr_Sgdinv)
//            due to the following nature A and B
//               A.  overlap matrix is positive define Hermite matrix
//               B.  inverse matrix of Hermite matrix is also Hermite matrix
//
// Caution  : 1. In this version, eigenvalue check is not supported.
//               so,naux must be set equall to naxucut.
//            2. This class is vailable only for real symmetric
//               positive define matrix
//            3. This class require squre full-matrix on memory
//               (half for MCD and the other for File blocking)
//
// Input    : overlap matrix (Hermite matrix)
//                fl_Mtr_Sab
//                fl_Mtr_Sgd
// Output   : inverse matrix (Hermite matrix)
//            (indexing : pq-through index,col oriented)
//            (format: index list + element,unpacked)
//                fl_Mtr_Sabinv
//                fl_Mtr_Sgdinv
//
// 暫定的に対角化ルーチンを利用して逆行列を作るようにした。


class DfInvMatrix : public DfObject {
public:
    DfInvMatrix(TlSerializeData* pPdfParam);
    virtual ~DfInvMatrix();

public:
    virtual void DfInvMain();

protected:
    template<class SymmetricMatrixType>
    void exec();
    
protected:
    bool m_bIsXcFitting;
};


template<class SymmetricMatrixType>
void DfInvMatrix::exec()
{
    // for Sab (cd basis)
    this->logger(" construct Sabinv with three full matrix in term of diagonalization");
    {
        SymmetricMatrixType Sab = DfObject::getSabMatrix<SymmetricMatrixType>();
        Sab.inverse();
        DfObject::saveSabInvMatrix(Sab);
    }
    
    // for Sgd (xc basis)
    if (this->m_bIsXcFitting == true) {
        this->logger(" construct Sgdinv with three full matrix in term of diagonalization");
        SymmetricMatrixType Sgd = DfObject::getSgdMatrix<SymmetricMatrixType>();
        Sgd.inverse();
        DfObject::saveSgdInvMatrix(Sgd);
    }
}


#endif // DFINVMATRIX_H

