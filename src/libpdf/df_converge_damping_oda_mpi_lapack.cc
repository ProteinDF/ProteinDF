#include "df_converge_damping_oda_mpi_lapack.h"

DfConverge_Damping_Oda_MpiLapack::DfConverge_Damping_Oda_MpiLapack(TlSerializeData* pPdfParam)
    : DfConverge_Damping_Oda_MpiTemplate<TlDenseSymmetricMatrix_Lapack, TlDenseVector_Lapack>(pPdfParam) {
}

DfConverge_Damping_Oda_MpiLapack::~DfConverge_Damping_Oda_MpiLapack() {}
