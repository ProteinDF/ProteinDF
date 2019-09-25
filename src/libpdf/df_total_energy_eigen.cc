#include "df_total_energy_eigen.h"

DfTotalEnergy_Eigen::DfTotalEnergy_Eigen(TlSerializeData* pPdfParam)
    : DfTotalEnergy_tmpl<TlDenseGeneralMatrix_Eigen,
                         TlDenseSymmetricMatrix_Eigen, TlDenseVector_Eigen,
                         DfOverlapX>(pPdfParam) {
    this->log_.info("run DfTotalEnergy_Eigen::DfTotalEnergy_Eigen()");
}

DfTotalEnergy_Eigen::~DfTotalEnergy_Eigen() {}
