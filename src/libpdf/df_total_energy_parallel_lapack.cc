#include "df_total_energy_parallel_lapack.h"

DfTotalEnergy_Parallel_Lapack::DfTotalEnergy_Parallel_Lapack(
    TlSerializeData* pPdfParam)
    : DfTotalEnergy_Parallel_tmpl<TlDenseGeneralMatrix_Lapack,
                                  TlDenseSymmetricMatrix_Lapack,
                                  TlDenseVector_Lapack, DfOverlapX>(pPdfParam) {
}

DfTotalEnergy_Parallel_Lapack::~DfTotalEnergy_Parallel_Lapack() {}
