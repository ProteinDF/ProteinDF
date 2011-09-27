#ifdef HAVE_CONFIG_H
#include "config.h"    // this file created by autotools
#endif // HAVE_CONFIG_H

#include "DfSummary_Parallel.h"
#include "TlDistributeMatrix.h"

DfSummary_Parallel::DfSummary_Parallel(TlSerializeData* pPdfParam)
    : DfSummary(pPdfParam)
{
}


DfSummary_Parallel::~DfSummary_Parallel()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    rComm.barrier();
}


void DfSummary_Parallel::logger(const std::string& str) const
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfObject::logger(str);
    }
}


void DfSummary_Parallel::exec()
{
#ifdef HAVE_SCALAPACK
    if (this->m_bUsingSCALAPACK == true) {
        this->exec_ScaLAPACK();
        return;
    }
#endif // HAVE_SCALAPACK

    this->exec_LAPACK();
}


void DfSummary_Parallel::exec_LAPACK()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    if (rComm.isMaster() == true) {
        DfSummary::exec();
    }
}


void DfSummary_Parallel::exec_ScaLAPACK()
{
    DfSummary::exec<TlDistributeMatrix>();
}
