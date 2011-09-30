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


// double DfTotalEnergy_Parallel::calculate_one_electron_part(const TlDistributeSymmetricMatrix& D)
// {
//     TlCommunicate& rComm = TlCommunicate::getInstance();

//     TlDistributeSymmetricMatrix Hpq = DfObject::getHpqMatrix<TlDistributeSymmetricMatrix>();
//     if (Hpq.getNumOfRows() != this->m_nNumOfAOs || Hpq.getNumOfCols() != this->m_nNumOfAOs) {
//         CnErr.abort("DfTotalenergy", "calculate_one_electron_part", "", "program error");
//     }

//     // add dummy charge
//     {
//         int number_dummyatom = 0;
//         if (rComm.isMaster() == true) {
//             Fl_Geometry geom(Fl_Geometry::getDefaultFileName());
//             number_dummyatom = geom.getDummyatom();
//         }
//         rComm.broadcast(number_dummyatom);

//         if (number_dummyatom != 0) {
//             const int chgextra_number = (*(this->pPdfParam_))["charge-extrapolate-number"].getInt();

//             TlDistributeSymmetricMatrix Hpq2 = DfObject::getHpq2Matrix<TlDistributeSymmetricMatrix>();

//             const int rotnum = std::min(this->m_nIteration, (chgextra_number +1));
//             Hpq += (rotnum * Hpq2);
//         }
//     }

//     double E_one_electron_part = 0.0;
//     const int dNumOfAOs = this->m_nNumOfAOs;
//     for (int i = 0; i < dNumOfAOs; ++i) {
//         // case: i != j
//         for (int j = 0; j < i; ++j) {
//             E_one_electron_part += 2.0 * D(i,j) * Hpq(i,j);
//         }
//         // case: i == j
//         E_one_electron_part += D(i,i) * Hpq(i,i);
//     }

//     return E_one_electron_part;
// }

