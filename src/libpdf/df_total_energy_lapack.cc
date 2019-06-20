#include "df_total_energy_lapack.h"

DfTotalEnergy_Lapack::DfTotalEnergy_Lapack(TlSerializeData* pPdfParam)
    : DfTotalEnergy_tmpl<TlDenseGeneralMatrix_Lapack,
                         TlDenseSymmetricMatrix_Lapack, TlDenseVector_Lapack,
                         DfOverlapX>(pPdfParam) {
    this->log_.info("run DfTotalEnergy_Lapack::DfTotalEnergy_Lapack()");
}

DfTotalEnergy_Lapack::~DfTotalEnergy_Lapack() {}
