#ifdef HAVE_CONFIG_H
#include "config.h"    // this file created by autotools
#endif // HAVE_CONFIG_H

#include <iostream>
#include <set>
#include <string>

#include "CnError.h"
#include "PdfUserInput.h"
#include "DfInitialguess.h"
#include "DfInitialGuessHarris.h"
#include "DfPopulation.h"

#include "DfOverlapX.h"
#include "TlSymmetricMatrix.h"
#include "TlSerializeData.h"

DfInitialGuessHarris::DfInitialGuessHarris(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam), debug_(false)
{
    std::string harrisDbFile = "harris.mpac";
    const char* pdfHome = std::getenv("PDF_HOME");
    if (pdfHome != NULL) {
        harrisDbFile = TlUtils::format("%s/data/harris.mpac", pdfHome);
    }

    TlMsgPack mpac;
    if (!mpac.load(harrisDbFile)) {
        this->log_.critical(TlUtils::format("cannot load %s.", harrisDbFile.c_str()));
        abort();
    }
    this->pdfParam_harrisDB_ = mpac.getSerializeData();
 
    if (this->isWorkOnDisk_ == true) {
        this->logger(" initial guess by harris functional is built on disk.\n");
        TlMatrix::useMemManager(true);
    }
}


DfInitialGuessHarris::~DfInitialGuessHarris()
{
}


void DfInitialGuessHarris::main()
{
    DfInitialGuessHarris::calcInitialDensityMatrix<TlMatrix, TlSymmetricMatrix, DfOverlapX, DfPopulation>();
}


