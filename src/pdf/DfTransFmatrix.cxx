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

#include "CnError.h"
#include "DfTransFmatrix.h"
#include "TlUtils.h"
#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"

DfTransFmatrix::DfTransFmatrix(TlSerializeData* pPdfParam, bool bExecDiis)
    : DfObject(pPdfParam), m_bExecDiis(bExecDiis)
{
}

DfTransFmatrix::~DfTransFmatrix()
{
}

void DfTransFmatrix::DfTrsFmatMain()
{
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
        CnErr.abort();
        break;
    }
}

void DfTransFmatrix::DfTrsFmatQclo(const std::string& fragname, int norbcut)
{
    this->m_nNumOfMOs = norbcut;

    switch (this->m_nMethodType) {
    case METHOD_RKS:
        this->main<TlMatrix, TlSymmetricMatrix>(RUN_RKS, fragname, true);
        break;

    case METHOD_UKS:
        this->main<TlMatrix, TlSymmetricMatrix>(RUN_UKS_ALPHA, fragname, true);
        this->main<TlMatrix, TlSymmetricMatrix>(RUN_UKS_BETA, fragname, true);
        break;

    case METHOD_ROKS:
        this->main<TlMatrix, TlSymmetricMatrix>(RUN_ROKS, fragname, true);
        break;

    default:
        CnErr.abort();
        break;
    }
}

