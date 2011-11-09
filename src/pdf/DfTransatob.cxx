#include "CnError.h"
#include "DfTransatob.h"
#include "TlUtils.h"
#include "TlMatrix.h"

DfTransatob::DfTransatob(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam)
{
//     this->scftype = flGbi["SCF"]["method"];
//     assert(this->scftype == "nsp" || this->scftype == "roks" || this->scftype == "sp");
}

DfTransatob::~DfTransatob()
{
}

void DfTransatob::DfTrsatobMain()
{
    switch (this->m_nMethodType) {
    case METHOD_RKS:
        this->main<TlMatrix>(RUN_RKS); // RKS
        break;

    case METHOD_UKS:
        this->main<TlMatrix>(RUN_UKS_ALPHA);  // UKS alpha spin
        this->main<TlMatrix>(RUN_UKS_BETA);  // UKS beta spin
        break;

    case METHOD_ROKS:
        this->main<TlMatrix>(RUN_ROKS);
        break;

    default:
        CnErr.abort();
        break;
    }
}

// for extended QCLO method
void DfTransatob::DfTrsatobQclo(const std::string& fragname, int norbcut)
{
    this->m_nNumOfMOs = norbcut;

    switch (this->m_nMethodType) {
    case METHOD_RKS:
        this->main<TlMatrix>(RUN_RKS, fragname, true); // RKS
        break;

    case METHOD_UKS:
        this->main<TlMatrix>(RUN_UKS_ALPHA, fragname, true);  // UKS alpha spin
        this->main<TlMatrix>(RUN_UKS_BETA, fragname, true);  // UKS beta spin
        break;

    case METHOD_ROKS:
        this->main<TlMatrix>(RUN_ROKS, fragname, true);
        break;

    default:
        CnErr.abort();
        break;
    }
}

