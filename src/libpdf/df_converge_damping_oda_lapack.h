#ifndef DF_CONVERGE_DAMPING_ODA_LAPACK_H
#define DF_CONVERGE_DAMPING_ODA_LAPACK_H

#include "df_converge_damping_oda_tmpl.h"
#include "tl_dense_symmetric_matrix_lapack.h"
#include "tl_dense_vector_lapack.h"

class DfConverge_Damping_Oda_Lapack
    : public DfConverge_Damping_Oda_tmpl<TlDenseSymmetricMatrix_Lapack, TlDenseVector_Lapack> {
   public:
    DfConverge_Damping_Oda_Lapack(TlSerializeData* pPdfParam);
    virtual ~DfConverge_Damping_Oda_Lapack();
};

#endif  // DF_CONVERGE_DAMPING_ODA_LAPACK_H
