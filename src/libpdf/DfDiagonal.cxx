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

#include "DfDiagonal.h"
#include "CnError.h"
#include "TlUtils.h"
#include "TlSymmetricMatrix.h"

DfDiagonal::DfDiagonal(TlSerializeData* pPdfParam) : DfObject(pPdfParam)
{
}


DfDiagonal::~DfDiagonal()
{
}

void DfDiagonal::DfDiagMain()
{
    // output informations
    switch (this->m_nMethodType) {
    case METHOD_RKS:
        this->main<TlMatrix, TlSymmetricMatrix>(RUN_RKS);
        break;

    case METHOD_UKS:
        this->main<TlMatrix, TlSymmetricMatrix>(RUN_UKS_ALPHA);
        this->main<TlMatrix, TlSymmetricMatrix>(RUN_UKS_BETA);
        break;

    case METHOD_ROKS:
        this->main<TlMatrix, TlSymmetricMatrix>(RUN_ROKS);
        break;

    default:
        CnErr.abort("DfDiagonal", "", "DfDiagMain", "the value of scftype is illegal");
        break;
    }
}


// for extended QCLO method
void DfDiagonal::DfDiagQclo(DfObject::RUN_TYPE runType, const std::string& fragname, int norbcut)
{
    this->m_nNumOfMOs = norbcut;
    this->main<TlMatrix, TlSymmetricMatrix>(runType, fragname, true);
}

