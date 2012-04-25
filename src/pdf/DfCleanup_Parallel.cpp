#include "TlCommunicate.h"
#include "DfCleanup_Parallel.h"

DfCleanup_Parallel::DfCleanup_Parallel(TlSerializeData* pPdfParam)
    : DfCleanup(pPdfParam) {
}


DfCleanup_Parallel::~DfCleanup_Parallel()
{
}


void DfCleanup_Parallel::cleanup()
{
    DfCleanup::cleanup();
}






