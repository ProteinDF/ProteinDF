#include "DfConverge.h"
#include "CnError.h"
#include "FileX.h"
#include "TlVector.h"
#include "TlSymmetricMatrix.h"
#include "TlUtils.h"

DfConverge::DfConverge(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam), m_nConvergeTarget(RHO_TILDE)
{
    const TlSerializeData& pdfParam = *pPdfParam;
    const std::string sDampingType  =
        TlUtils::toUpper(pdfParam["model"]["scf-acceleration/damping/damping-type"].getStr());
    if (sDampingType == "DENSITY") {
        this->m_nConvergeTarget = RHO_TILDE;
    } else if (sDampingType == "FOCK") {
        this->m_nConvergeTarget = KS_MATRIX;
    } else if (sDampingType == "DENSITY_MATRIX") {
        this->m_nConvergeTarget = DENSITY_MATRIX;
    } else {
        CnErr.abort("unknown damping type. stop.");
    }
}

DfConverge::~DfConverge()
{
}

void DfConverge::doConverge()
{
    switch (this->m_nConvergeTarget) {
    case RHO_TILDE:
        this->convergeRhoTilde();
        break;
    case KS_MATRIX:
        this->convergeKSMatrix();
        break;
    case DENSITY_MATRIX:
        this->convergePMatrix();
        break;
    default:
        CnErr.abort("Unknown converge target. stop.");
        break;
    }
}

