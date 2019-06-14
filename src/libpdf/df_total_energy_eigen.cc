#include "df_total_energy_eigen.h"

DfTotalEnergy_Eigen::DfTotalEnergy_Eigen(TlSerializeData* pPdfParam)
    : DfTotalEnergy_tmpl<
          TlDenseGeneralMatrix_Eigen, TlDenseSymmetricMatrix_Eigen,
          TlDenseVector_Eigen, DfOverlapX>(pPdfParam) {}

DfTotalEnergy_Eigen::~DfTotalEnergy_Eigen() {}
