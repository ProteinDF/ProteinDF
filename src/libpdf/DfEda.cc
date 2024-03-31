#include "DfEda.h"

#include "DfGridFreeXC.h"

DfEda::DfEda(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam) {
}

DfEda::~DfEda() {}

void DfEda::calc() {
    if (this->isDFT_ == true) {
        switch (this->XC_engine_) {
            case XC_ENGINE_GRID:
                this->calcGridFreeXC();
                break;

            case XC_ENGINE_GRIDFREE:
            case XC_ENGINE_GRIDFREE_CD:
                // DO NOTHING
                break;

            default:
                this->log_.critical("program error.");
                break;
        }
    }
}

void DfEda::calcGridFreeXC() {
    DfGridFreeXC* pDfGridFreeXC = this->getDfGridFreeXcObject();

    this->log_.info("prepare gridfree XC term");
    pDfGridFreeXC->preprocessBeforeSCF();

    this->log_.info("calc gridfree XC term");
    pDfGridFreeXC->buildFxc();

    delete pDfGridFreeXC;
    pDfGridFreeXC = NULL;
}

DfGridFreeXC* DfEda::getDfGridFreeXcObject() {
    DfGridFreeXC* pDfGridFreeXC = new DfGridFreeXC(this->pPdfParam_);
    return pDfGridFreeXC;
}
