#include "df_converge_damping_diis_mpi_lapack.h"

DfConverge_Damping_Diis_Mpi_Lapack::DfConverge_Damping_Diis_Mpi_Lapack(TlSerializeData* pPdfParam)
    : DfConverge_Damping_Diis_Mpi_Template<TlDenseGeneralMatrix_Lapack, TlDenseSymmetricMatrix_Lapack, TlDenseVector_Lapack>(pPdfParam) {
}

DfConverge_Damping_Diis_Mpi_Lapack::~DfConverge_Damping_Diis_Mpi_Lapack() {
}
