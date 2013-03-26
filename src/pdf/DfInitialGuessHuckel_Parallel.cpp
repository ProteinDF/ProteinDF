#include <iostream>
#include <limits>
#include "DfInitialGuessHuckel_Parallel.h"
#include "TlDistributeMatrix.h"
#include "TlDistributeSymmetricMatrix.h"

DfInitialGuessHuckel_Parallel::DfInitialGuessHuckel_Parallel(TlSerializeData* pPdfParam)
    : DfInitialGuessHuckel(pPdfParam)
{
}

DfInitialGuessHuckel_Parallel::~DfInitialGuessHuckel_Parallel()
{
}

void DfInitialGuessHuckel_Parallel::createGuess()
{
    DfInitialGuessHuckel::createGuess<TlDistributeSymmetricMatrix, TlDistributeMatrix>();
}

