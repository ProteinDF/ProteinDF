#ifndef DF_CONVERGE_DAMPING_ANDERSON_LAPACK_H
#define DF_CONVERGE_DAMPING_ANDERSON_LAPACK_H

#include "df_converge_damping_anderson_template.h"
#include "tl_dense_general_matrix_lapack.h"
#include "tl_dense_symmetric_matrix_lapack.h"
#include "tl_dense_vector_lapack.h"

class DfConverge_Damping_Anderson_Lapack : public DfConverge_Damping_Anderson_Template<TlDenseGeneralMatrix_Lapack, TlDenseSymmetricMatrix_Lapack, TlDenseVector_Lapack> {
public:
    DfConverge_Damping_Anderson_Lapack(TlSerializeData* pPdfParam);
    virtual ~DfConverge_Damping_Anderson_Lapack();
};

#endif  // DF_CONVERGE_DAMPING_ANDERSON_LAPACK_H
