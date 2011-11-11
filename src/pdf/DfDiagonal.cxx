#ifdef HAVE_CONFIG_H
#include "config.h"    // this file created by autotools
#endif // HAVE_CONFIG_H

#include "DfDiagonal.h"
#include "CnError.h"
#include "TlUtils.h"
#include "TlSymmetricMatrix.h"

DfDiagonal::DfDiagonal(TlSerializeData* pPdfParam) : DfObject(pPdfParam)
{
}


DfDiagonal::~DfDiagonal()
{
}

void DfDiagonal::DfDiagMain()
{
    // output informations
    switch (this->m_nMethodType) {
    case METHOD_RKS:
        this->main<TlMatrix, TlSymmetricMatrix>(RUN_RKS);
        break;

    case METHOD_UKS:
        this->main<TlMatrix, TlSymmetricMatrix>(RUN_UKS_ALPHA);
        this->main<TlMatrix, TlSymmetricMatrix>(RUN_UKS_BETA);
        break;

    case METHOD_ROKS:
        this->main<TlMatrix, TlSymmetricMatrix>(RUN_ROKS);
        break;

    default:
        CnErr.abort("DfDiagonal", "", "DfDiagMain", "the value of scftype is illegal");
        break;
    }
}


// for extended QCLO method
void DfDiagonal::DfDiagQclo(DfObject::RUN_TYPE runType, const std::string& fragname, int norbcut)
{
    this->m_nNumOfMOs = norbcut;
    this->main<TlMatrix, TlSymmetricMatrix>(runType, fragname, true);
}

