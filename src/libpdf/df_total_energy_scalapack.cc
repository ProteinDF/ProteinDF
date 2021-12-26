#include "df_total_energy_scalapack.h"

#include "TlCommunicate.h"

DfTotalEnergy_Scalapack::DfTotalEnergy_Scalapack(TlSerializeData* pPdfParam)
    : DfTotalEnergyObject(pPdfParam) {
    this->log_.info("run DfTotalEnergy_Scalapack::DfTotalEnergy_Scalapack()");
}

DfTotalEnergy_Scalapack::~DfTotalEnergy_Scalapack() {}

void DfTotalEnergy_Scalapack::output_stdout(const std::string& str) const {
    TlCommunicate& rComm = TlCommunicate::getInstance();
    if (rComm.isMaster()) {
        DfTotalEnergyObject::output_stdout(str);
    }
}
