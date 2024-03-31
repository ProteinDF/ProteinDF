#ifndef DF_CONVERGE_DAMPING_MPI_LAPACK_H
#define DF_CONVERGE_DAMPING_MPI_LAPACK_H

#include "df_converge_damping_mpi_template.h"
#include "tl_dense_symmetric_matrix_lapack.h"
#include "tl_dense_vector_lapack.h"

class DfConverge_Damping_Mpi_Lapack : public DfConverge_Damping_Mpi_Template<TlDenseSymmetricMatrix_Lapack, TlDenseVector_Lapack> {
public:
    DfConverge_Damping_Mpi_Lapack(TlSerializeData* pPdfParam);
    virtual ~DfConverge_Damping_Mpi_Lapack();
};

#endif  // DF_CONVERGE_DAMPING_MPI_LAPACK_H
