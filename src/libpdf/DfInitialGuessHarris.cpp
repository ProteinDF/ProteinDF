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
#include "config.h"  // this file created by autotools
#endif               // HAVE_CONFIG_H

#include <iostream>
#include <set>
#include <string>

#include "CnError.h"
#include "DfInitialGuessHarris.h"
#include "DfPopulation.h"
#include "PdfUserInput.h"

#include "DfOverlapX.h"
#include "TlSerializeData.h"
#include "tl_dense_symmetric_matrix_lapack.h"

DfInitialGuessHarris::DfInitialGuessHarris(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam), debug_(false) {
    // if (this->isWorkOnDisk_ == true) {
    //   this->logger(" initial guess by harris functional is built on
    //   disk.\n"); TlDenseGeneralMatrix_Lapack::useMemManager(true);
    // }
}

DfInitialGuessHarris::~DfInitialGuessHarris() {}

void DfInitialGuessHarris::main() {
    this->loadHarrisDB();

    DfInitialGuessHarris::calcInitialDensityMatrix<
        TlDenseGeneralMatrix_Lapack, TlDenseSymmetricMatrix_Lapack, DfOverlapX,
        DfPopulation>();
}

void DfInitialGuessHarris::loadHarrisDB() {
    std::string harrisDbFile = "harris.mpac";
    const char* pdfHome = std::getenv("PDF_HOME");
    if (pdfHome != NULL) {
        harrisDbFile = TlUtils::format("%s/data/harris.mpac", pdfHome);
    }

    TlMsgPack mpac;
    if (!mpac.load(harrisDbFile)) {
        this->log_.critical(
            TlUtils::format("cannot load %s.", harrisDbFile.c_str()));
        abort();
    }
    this->pdfParam_harrisDB_ = mpac.getSerializeData();
}
