#ifdef HAVE_CONFIG_H
#include "config.h"    // this file created by autotools
#endif // HAVE_CONFIG_H

#include "DfInvMatrix_Parallel.h"
#include "TlCommunicate.h"
#include "TlSymmetricMatrix.h"
#include "TlDistributeSymmetricMatrix.h"

DfInvMatrix_Parallel::DfInvMatrix_Parallel(TlSerializeData* pPdfParam)
    : DfInvMatrix(pPdfParam)
{
}


DfInvMatrix_Parallel::~DfInvMatrix_Parallel()
{
}


void DfInvMatrix_Parallel::logger(const std::string& str) const
{
    TlCommunicate& rComm = TlCommunicate::getInstance();

    if (rComm.isMaster() == true) {
        DfInvMatrix::logger(str);
    }
}


void DfInvMatrix_Parallel::DfInvMain()
{
#ifdef HAVE_SCALAPACK
    if (this->m_bUsingSCALAPACK == true) {
        this->exec<TlDistributeSymmetricMatrix>();
        return;
    }
#endif // HAVE_SCALAPACK

    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster() == true) {
        DfInvMatrix::DfInvMain();
    }
    rComm.barrier();
}


// void DfInvMatrix_Parallel::inverseMatrix(const std::string& sInFile, const std::string& sOutFile)
// {
//     TlCommunicate& rComm = TlCommunicate::getInstance();

// #ifdef HAVE_SCALAPACK
//     if (this->m_bUsingSCALAPACK == true) {
//         this->inverseMatrix_ScaLAPACK(sInFile, sOutFile);
//     } else {
//         this->inverseMatrix_LAPACK(sInFile, sOutFile);
//     }
// #else
//     this->inverseMatrix_LAPACK(sInFile, sOutFile);
// #endif // HAVE_SCALAPACK

//     rComm.barrier();
// }


// void DfInvMatrix_Parallel::inverseMatrix_LAPACK(const std::string& sInFile,
//                                                 const std::string& sOutFile)
// {
//     this->logger(" Inverse matrix is operated using LAPACK.\n");
//     TlCommunicate& rComm = TlCommunicate::getInstance();
//     if (rComm.isMaster() == true) {
//         DfInvMatrix::inverseMatrix_tmpl<TlSymmetricMatrix>(sInFile, sOutFile);
//     }
// }


// void DfInvMatrix_Parallel::inverseMatrix_ScaLAPACK(const std::string& sInFile,
//                                                    const std::string& sOutFile)
// {
//     this->logger(" Inverse matrix is operated using ScaLAPACK.\n");
//     DfInvMatrix::inverseMatrix_tmpl<TlDistributeSymmetricMatrix>(sInFile, sOutFile);
// }


