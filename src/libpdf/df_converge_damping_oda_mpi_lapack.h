#ifndef DF_CONVERGE_DAMPING_ODA_MPI_LAPACK_H
#define DF_CONVERGE_DAMPING_ODA_MPI_LAPACK_H

#include "df_converge_damping_oda_mpi_template.h"
#include "tl_dense_symmetric_matrix_lapack.h"
#include "tl_dense_vector_lapack.h"

class DfConverge_Damping_Oda_MpiLapack : public DfConverge_Damping_Oda_MpiTemplate<TlDenseSymmetricMatrix_Lapack, TlDenseVector_Lapack> {
public:
    DfConverge_Damping_Oda_MpiLapack(TlSerializeData* pPdfParam);
    virtual ~DfConverge_Damping_Oda_MpiLapack();
};

#endif  // DF_CONVERGE_DAMPING_ODA_MPI_LAPACK_H
