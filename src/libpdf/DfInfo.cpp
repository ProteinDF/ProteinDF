#include "DfInfo.h"
#include <iostream>

DfInfo::DfInfo(TlSerializeData* pPdfParam) : DfObject(pPdfParam) {}

DfInfo::~DfInfo() {}

int DfInfo::getNumOfElectrons() const {
    int numOfElectrons = 0;
    switch (this->m_nMethodType) {
        case METHOD_RKS:
            numOfElectrons = this->m_nNumOfElectrons;
            break;

        case METHOD_UKS:
            numOfElectrons = this->m_nNumOfAlphaElectrons;
            numOfElectrons += this->m_nNumOfBetaElectrons;
            break;

        case METHOD_ROKS:
            numOfElectrons = this->numOfClosedShellElectrons_;
            numOfElectrons += this->numOfOpenShellElectrons_;
            break;

        default:
            std::cerr << "unknown method." << std::endl;
            break;
    }

    return numOfElectrons;
}
