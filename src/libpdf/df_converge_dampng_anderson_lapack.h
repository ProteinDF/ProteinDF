#ifndef DF_CONVERGE_ANDERSON_LAPACK_H
#define DF_CONVERGE_ANDERSON_LAPACK_H

#include "df_converge_anderson_template.h"
#include "tl_dense_general_matrix_lapack.h"
#include "tl_dense_symmetric_matrix_lapack.h"
#include "tl_dense_vector_lapack.h"

class DfConverge_Anderson_Lapack : public DfConverge_AndersonTemplate<TlDenseGeneralMatrix_Lapack, TlDenseSymmetricMatrix_Lapack, TlDenseVector_Lapack> {
public:
    DfConverge_Anderson_Lapack(TlSerializeData* pPdfParam);
    virtual ~DfConverge_Anderson_Lapack();
};

#endif  // DF_CONVERGE_ANDERSON_LAPACK_H
