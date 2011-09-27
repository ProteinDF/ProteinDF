#include "DfJMatrix.h"
#include "DfEriX.h"
#include "TlSymmetricMatrix.h"
#include "TlUtils.h"

DfJMatrix::DfJMatrix(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam) {
}


DfJMatrix::~DfJMatrix()
{
}


void DfJMatrix::buildJMatrix()
{
    const index_type numOfAOs = this->m_nNumOfAOs;
    const int iteration = this->m_nIteration;
    const RUN_TYPE runType = RUN_RKS;
    
    const TlSymmetricMatrix P = this->getDiffDensityMatrix<TlSymmetricMatrix>(runType, iteration);

    DfEriX dfEri(this->pPdfParam_);
    TlSymmetricMatrix J(numOfAOs);
    dfEri.getJpq(P, &J);

    if (iteration >= 2) {
        const TlSymmetricMatrix prevJ = this->getJMatrix<TlSymmetricMatrix>(iteration -1);
        J += prevJ;
    }
    
    this->saveJMatrix(iteration, J);
}



