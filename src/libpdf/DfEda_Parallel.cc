#include "DfEda_Parallel.h"

#include "DfGridFreeXC_Parallel.h"

DfEda_Parallel::DfEda_Parallel(TlSerializeData* pPdfParam)
    : DfEda(pPdfParam) {
}

DfEda_Parallel::~DfEda_Parallel() {
}

DfGridFreeXC* DfEda_Parallel::getDfGridFreeXcObject() {
    return new DfGridFreeXC_Parallel(this->pPdfParam_);
}
