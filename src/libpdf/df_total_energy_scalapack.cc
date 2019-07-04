#include "df_total_energy_scalapack.h"

DfTotalEnergy_Scalapack::DfTotalEnergy_Scalapack(TlSerializeData* pPdfParam)
    : DfTotalEnergy_tmpl<TlDenseGeneralMatrix_Scalapack,
                         TlDenseSymmetricMatrix_Scalapack,
                         TlDenseVector_Scalapack, DfOverlapX_Parallel>(
          pPdfParam) {
    this->log_.info("run DfTotalEnergy_Scalapack::DfTotalEnergy_Scalapack()");
          }

DfTotalEnergy_Scalapack::~DfTotalEnergy_Scalapack() {}
