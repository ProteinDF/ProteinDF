#include <cassert>
#include <fstream>
#include <string>
#include <cmath>

#include "DfXMatrix.h"
#include "TlVector.h"
#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"

DfXMatrix::DfXMatrix(TlSerializeData* pPdfParam) : DfObject(pPdfParam)
{
    assert(pPdfParam != NULL);
    const TlSerializeData& pdfParam = *pPdfParam;
    this->threshold_trancation = pdfParam["orbital-independence-threshold"].getDouble();

    this->debug_save_mat_ = pdfParam["debug/DfXMatrix/save_mat"].getBoolean();
    this->debug_check_X_ = pdfParam["debug/DfXMatrix/check_X"].getBoolean();
}

DfXMatrix::~DfXMatrix()
{
}

void DfXMatrix::buildX()
{
    TlSymmetricMatrix S = this->getSpqMatrix<TlSymmetricMatrix>();
    TlMatrix X;
    TlMatrix Xinv;
    
    this->canonicalOrthogonalize(S, &X, &Xinv);

    DfObject::saveXMatrix(X);
    DfObject::saveXInvMatrix(Xinv);
    (*(this->pPdfParam_))["num_of_MOs"] = X.getNumOfCols();
}


