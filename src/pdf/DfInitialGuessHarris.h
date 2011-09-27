#ifndef DFINITIALGUESSHARRIS_H
#define DFINITIALGUESSHARRIS_H

#include <string>
#include "DfObject.h"
#include "TlVector.h"
#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"

/// Harrisの汎関数による初期値を作成する
class DfInitialGuessHarris : public DfObject {
private:
    enum ScfType {
        RKS,
        UKS,
        ROKS
    };

public:
    DfInitialGuessHarris(TlSerializeData* pPdfParam);
    virtual ~DfInitialGuessHarris();

    void main();
    
private:
    TlSerializeData pdfParam_harrisDB_;
};

#endif // DFINITIALGUESSHarris_H
