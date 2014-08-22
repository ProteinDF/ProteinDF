// Copyright (C) 2002-2014 The ProteinDF project
// see also AUTHORS and README.
// 
// This file is part of ProteinDF.
// 
// ProteinDF is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// ProteinDF is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

#ifdef HAVE_CONFIG_H
#include "config.h"    // this file created by autotools
#endif // HAVE_CONFIG_H

#include <iostream>
#include <set>
#include <string>

#include "CnError.h"
#include "PdfUserInput.h"
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


