#ifndef DFINITIALGUESSHUCKEL_PARALLEL_H
#define DFINITIALGUESSHUCKEL_PARALLEL_H

#include <string>
#include "DfInitialGuessHuckel.h"

class DfInitialGuessHuckel_Parallel : public DfInitialGuessHuckel {
public:
    DfInitialGuessHuckel_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfInitialGuessHuckel_Parallel();

public:
    virtual void createGuess();
    
};

#endif // DFINITIALGUESSHUCKEL_PARALLEL_H
