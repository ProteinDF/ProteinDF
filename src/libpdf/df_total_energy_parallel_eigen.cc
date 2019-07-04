#include "df_total_energy_parallel_eigen.h"

DfTotalEnergy_Parallel_Eigen::DfTotalEnergy_Parallel_Eigen(
    TlSerializeData* pPdfParam)
    : DfTotalEnergy_Parallel_tmpl<TlDenseGeneralMatrix_Eigen,
                                  TlDenseSymmetricMatrix_Eigen,
                                  TlDenseVector_Eigen, DfOverlapX>(pPdfParam) {
}

DfTotalEnergy_Parallel_Eigen::~DfTotalEnergy_Parallel_Eigen() {}
