#ifndef DF_TOTAL_ENERGY_PARALLEL_LAPACK_H
#define DF_TOTAL_ENERGY_PARALLEL_LAPACK_H

#include "df_total_energy_parallel_tmpl.h"

#include "DfOverlapX.h"
#include "tl_dense_general_matrix_lapack.h"
#include "tl_dense_symmetric_matrix_lapack.h"
#include "tl_dense_vector_lapack.h"

class DfTotalEnergy_Parallel_Lapack
    : public DfTotalEnergy_Parallel_tmpl<TlDenseGeneralMatrix_Lapack,
                                         TlDenseSymmetricMatrix_Lapack,
                                         TlDenseVector_Lapack, DfOverlapX> {
   public:
    DfTotalEnergy_Parallel_Lapack(TlSerializeData* pPdfParam);
    virtual ~DfTotalEnergy_Parallel_Lapack();
};

#endif  // DF_TOTAL_ENERGY_PARALLEL_LAPACK_H