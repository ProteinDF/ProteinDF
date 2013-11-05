#include "DfInvMatrix.h"
#include "TlSymmetricMatrix.h"

DfInvMatrix::DfInvMatrix(TlSerializeData* pPdfParam) : DfObject(pPdfParam)
{
    const TlSerializeData& pdfParam = *pPdfParam;
    const std::string sXcPotential = TlUtils::toUpper(pdfParam["xc_functional"].getStr());
    const char nLastChar = sXcPotential[sXcPotential.length() -1];
    this->m_bIsXcFitting = (nLastChar == '~') ? true : false;
}


DfInvMatrix::~DfInvMatrix()
{
}


void DfInvMatrix::DfInvMain()
{
    this->exec<TlSymmetricMatrix>();
}

