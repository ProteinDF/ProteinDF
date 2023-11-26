#include "df_converge_damping_mpi_lapack.h"

DfConverge_Damping_Mpi_Lapack::DfConverge_Damping_Mpi_Lapack(TlSerializeData* pPdfParam)
    : DfConverge_Damping_Mpi_Template<TlDenseSymmetricMatrix_Lapack, TlDenseVector_Lapack>(pPdfParam) {
}

DfConverge_Damping_Mpi_Lapack::~DfConverge_Damping_Mpi_Lapack() {}
