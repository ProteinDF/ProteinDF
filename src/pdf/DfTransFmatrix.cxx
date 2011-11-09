#include "CnError.h"
#include "DfTransFmatrix.h"
#include "TlUtils.h"
#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"

DfTransFmatrix::DfTransFmatrix(TlSerializeData* pPdfParam, bool bExecDiis)
    : DfObject(pPdfParam), m_bExecDiis(bExecDiis)
{
}

DfTransFmatrix::~DfTransFmatrix()
{
}

void DfTransFmatrix::DfTrsFmatMain()
{
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
        CnErr.abort();
        break;
    }
}

void DfTransFmatrix::DfTrsFmatQclo(const std::string& fragname, int norbcut)
{
    this->m_nNumOfMOs = norbcut;

    switch (this->m_nMethodType) {
    case METHOD_RKS:
        this->main<TlMatrix, TlSymmetricMatrix>(RUN_RKS, fragname, true);
        break;

    case METHOD_UKS:
        this->main<TlMatrix, TlSymmetricMatrix>(RUN_UKS_ALPHA, fragname, true);
        this->main<TlMatrix, TlSymmetricMatrix>(RUN_UKS_BETA, fragname, true);
        break;

    case METHOD_ROKS:
        this->main<TlMatrix, TlSymmetricMatrix>(RUN_ROKS, fragname, true);
        break;

    default:
        CnErr.abort();
        break;
    }
}

