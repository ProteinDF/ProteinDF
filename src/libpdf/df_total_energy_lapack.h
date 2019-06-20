#ifndef DF_TOTAL_ENERGY_LAPACK_H
#define DF_TOTAL_ENERGY_LAPACK_H

#include "DfOverlapX.h"
#include "df_total_energy_tmpl.h"
#include "tl_dense_general_matrix_lapack.h"
#include "tl_dense_symmetric_matrix_lapack.h"
#include "tl_dense_vector_lapack.h"

class DfTotalEnergy_Lapack
    : public DfTotalEnergy_tmpl<TlDenseGeneralMatrix_Lapack,
                                TlDenseSymmetricMatrix_Lapack,
                                TlDenseVector_Lapack, DfOverlapX> {
   public:
    DfTotalEnergy_Lapack(TlSerializeData* pPdfParam);
    virtual ~DfTotalEnergy_Lapack();
};

#endif  // DF_TOTAL_ENERGY_LAPACK_H
