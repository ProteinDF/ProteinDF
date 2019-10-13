#ifndef DF_TOTAL_ENERGY_SCALAPACK_H
#define DF_TOTAL_ENERGY_SCALAPACK_H

#include "DfOverlapX_Parallel.h"
#include "df_total_energy_tmpl.h"
#include "tl_dense_general_matrix_scalapack.h"
#include "tl_dense_symmetric_matrix_scalapack.h"
#include "tl_dense_vector_scalapack.h"

class DfTotalEnergy_Scalapack
    : public DfTotalEnergy_tmpl<TlDenseGeneralMatrix_Scalapack,
                                TlDenseSymmetricMatrix_Scalapack,
                                TlDenseVector_Scalapack, DfOverlapX_Parallel> {
   public:
    DfTotalEnergy_Scalapack(TlSerializeData* pPdfParam);
    virtual ~DfTotalEnergy_Scalapack();
};

#endif  // DF_TOTAL_ENERGY_SCALAPACK_H
