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

#include "df_converge.h"

#include "CnError.h"
#include "TlUtils.h"

DfConverge::DfConverge(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam), convergeTarget_(RHO_TILDE) {
    const TlSerializeData& pdfParam = *pPdfParam;
    const std::string dampingType = TlUtils::toUpper(pdfParam["scf_acceleration/damping/damping_type"].getStr());
    if (dampingType == "DENSITY") {
        this->convergeTarget_ = RHO_TILDE;
    } else if (dampingType == "FOCK") {
        this->convergeTarget_ = KS_MATRIX;
    } else if (dampingType == "DENSITY_MATRIX") {
        this->convergeTarget_ = DENSITY_MATRIX;
    } else {
        CnErr.abort("unknown damping type. stop.");
    }
}

DfConverge::~DfConverge() {}

void DfConverge::doConverge() {
    switch (this->convergeTarget_) {
        case RHO_TILDE:
            this->convergeRhoTilde();
            break;
        case KS_MATRIX:
            this->convergeKSMatrix();
            break;
        case DENSITY_MATRIX:
            this->convergePMatrix();
            break;
        default:
            CnErr.abort("Unknown converge target. stop.");
            break;
    }
}
