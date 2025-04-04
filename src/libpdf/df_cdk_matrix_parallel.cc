#ifdef HAVE_CONFIG_H
#include "config.h"  // this file created by autotools
#endif               // HAVE_CONFIG_H

#include "DfTaskCtrl_Parallel.h"
#include "TlCommunicate.h"
#include "df_cdk_matrix_parallel.h"

#ifdef HAVE_SCALAPACK
#include "tl_dense_general_matrix_scalapack.h"
#include "tl_dense_symmetric_matrix_scalapack.h"
#include "tl_dense_vector_scalapack.h"
#endif  // HAVE_SCALAPACK

// ----------------------------------------------------------------------------
// construct & destruct
// ----------------------------------------------------------------------------
DfCdkMatrix_Parallel::DfCdkMatrix_Parallel(TlSerializeData* pPdfParam)
    : DfCdkMatrix(pPdfParam) {}

DfCdkMatrix_Parallel::~DfCdkMatrix_Parallel() {}

// ----------------------------------------------------------------------------
// public
// ----------------------------------------------------------------------------
void DfCdkMatrix_Parallel::getK() {
    switch (this->linearAlgebraPackage_) {
#ifdef HAVE_SCALAPACK
        case LAP_SCALAPACK: {
            this->log_.info("Linear Algebra Package: ScaLapack");
            this->getK_runType<TlDenseSymmetricMatrix_Scalapack,
                               TlDenseVector_Scalapack,
                               TlDenseGeneralMatrix_Scalapack,
                               TlDenseSymmetricMatrix_Scalapack,
                               TlDenseSymmetricMatrix_Scalapack>();

        } break;
#endif  // HAVE_SCALAPACK

        default:
            this->log_.info("parallel CD-K matrix calculation");
            DfCdkMatrix::getK();
            break;
    }
}

// -----------------------------------------------------------------------------
// task control
// -----------------------------------------------------------------------------
DfTaskCtrl* DfCdkMatrix_Parallel::getDfTaskCtrlObject() const {
    this->log_.info("task control: parallel");
    DfTaskCtrl* pDfTaskCtrl = new DfTaskCtrl_Parallel(this->pPdfParam_);
    return pDfTaskCtrl;
}

// void DfCdkMatrix_Parallel::finalize(TlDenseSymmetricMatrixObject* pK) {
//     TlCommunicate& rComm = TlCommunicate::getInstance();
//     rComm.allReduce_SUM(*pK);
// }
