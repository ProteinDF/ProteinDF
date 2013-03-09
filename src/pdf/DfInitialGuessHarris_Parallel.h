#ifndef DFINITIALGUESSHARRIS_PARALLEL_H
#define DFINITIALGUESSHARRIS_PARALLEL_H

#include "DfInitialGuessHarris.h"

/// Harrisの汎関数による初期値を作成する
class DfInitialGuessHarris_Parallel : public DfInitialGuessHarris {
public:
    DfInitialGuessHarris_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfInitialGuessHarris_Parallel();

    virtual void main();
};

#endif // DFINITIALGUESSHARRIS_PARALLEL_H
