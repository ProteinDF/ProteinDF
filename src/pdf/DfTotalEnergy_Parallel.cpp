#ifdef HAVE_CONFIG_H
#include "config.h"    // this file created by autotools
#endif // HAVE_CONFIG_H

#include "DfTotalEnergy_Parallel.h"
#include "DfEri_Parallel.h"
#include "DfOverlap_Parallel.h"
#include "DfXCFunctional_Parallel.h"
#include "Fl_Geometry.h"
#include "TlCommunicate.h"
#include "CnError.h"

DfTotalEnergy_Parallel::DfTotalEnergy_Parallel(TlSerializeData* pPdfParam)
    : DfTotalEnergy(pPdfParam)
{
    this->m_bUseDistributeMatrix = false;
}


DfTotalEnergy_Parallel::~DfTotalEnergy_Parallel()
{
}


void DfTotalEnergy_Parallel::logger(const std::string& str) const
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    if (rComm.isMaster() == true) {
        DfTotalEnergy::logger(str);
    }
}


void DfTotalEnergy_Parallel::output()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfTotalEnergy::output();
    }
}


void DfTotalEnergy_Parallel::exec()
{
#ifdef HAVE_SCALAPACK
    if (this->m_bUsingSCALAPACK == true) {
        this->exec_ScaLAPACK();
        return;
    }
#endif // HAVE_SCALAPACK

    this->exec_LAPACK();
}


void DfTotalEnergy_Parallel::exec_LAPACK()
{
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfTotalEnergy::exec();
    }
    rComm.barrier();
}


void DfTotalEnergy_Parallel::exec_ScaLAPACK()
{
    DfTotalEnergy::exec_template<DfOverlap_Parallel, DfEri_Parallel, TlDistributeSymmetricMatrix, TlDistributeVector>();
}


// total energy including dummy charge
void DfTotalEnergy_Parallel::calculate_real_energy()
{
#ifdef HAVE_SCALAPACK
    if (this->m_bUsingSCALAPACK == true) {
        // ScaLAPACK
        this->calcRealEnergy<TlDistributeSymmetricMatrix>();
        return;
    }
#endif // HAVE_SCALAPACK

    // LAPACK
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfTotalEnergy::calculate_real_energy();
    }
}


