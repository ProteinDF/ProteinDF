#ifndef DFPOPULATION_PARALLEL_H
#define DFPOPULATION_PARALLEL_H

#include "DfPopulation.h"

class DfPopulation_Parallel : public DfPopulation {
public:
    DfPopulation_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfPopulation_Parallel();

protected:
    virtual void calcPop(const int iteration);
};

#endif // DFPOPULATION_PARALLEL_H
