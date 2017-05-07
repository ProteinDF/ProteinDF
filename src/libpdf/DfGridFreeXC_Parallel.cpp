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

#include "DfGridFreeXC_Parallel.h"
#include "DfOverlapX_Parallel.h"
#include "DfXMatrix_Parallel.h"
#include "DfCD_Parallel.h"

DfGridFreeXC_Parallel::DfGridFreeXC_Parallel(TlSerializeData* pPdfParam)
    : DfGridFreeXC(pPdfParam) {
    this->log_.info("DfGridFreeXC_Parallel::DfGridFreeXC_Parallel()");
}

DfGridFreeXC_Parallel::~DfGridFreeXC_Parallel()
{
}

DfOverlapX* DfGridFreeXC_Parallel::getDfOverlapObject()
{
    DfOverlapX* pDfOverlapX = new DfOverlapX_Parallel(this->pPdfParam_);
    return pDfOverlapX;
}

DfXMatrix* DfGridFreeXC_Parallel::getDfXMatrixObject()
{
    DfXMatrix* pDfXMatrix = new DfXMatrix_Parallel(this->pPdfParam_);
    return pDfXMatrix;
}

// before SCF ==================================================================
void DfGridFreeXC_Parallel::preprocessBeforeSCF()
{
#ifdef HAVE_SCALAPACK
    if (this->m_bUsingSCALAPACK == true) {
        this->preprocessBeforeSCF_ScaLAPACK();
    } else {
        this->preprocessBeforeSCF_LAPACK();
    }
#else
    {
        this->preprocessBeforeSCF_LAPACK();
    }
#endif // HAVE_SCALAPACK
}

void DfGridFreeXC_Parallel::preprocessBeforeSCF_LAPACK()
{
    this->preprocessBeforeSCF_templ<DfOverlapX_Parallel,
                                    DfXMatrix_Parallel,
                                    TlSymmetricMatrix, TlMatrix>();
}

void DfGridFreeXC_Parallel::preprocessBeforeSCF_ScaLAPACK()
{
    this->preprocessBeforeSCF_templ<DfOverlapX_Parallel,
                                    DfXMatrix_Parallel,
                                    TlDistributeSymmetricMatrix,
                                    TlDistributeMatrix>();
}


// in SCF ======================================================================
void DfGridFreeXC_Parallel::buildFxc_LDA()
{
    this->log_.info("DfGridFreeXC_Parallel::buildFxc_LDA()");
#ifdef HAVE_SCALAPACK
    if (this->m_bUsingSCALAPACK == true) {
        this->buildFxc_LDA_ScaLAPACK();
    } else {
        this->buildFxc_LDA_LAPACK();
    }
#else
    {
        this->buildFxc_LDA_LAPACK();
    }
#endif // HAVE_SCALAPACK
}

void DfGridFreeXC_Parallel::buildFxc_LDA_LAPACK()
{
    this->log_.info("DfGridFreeXC_Parallel::buildFxc_LDA_LAPACK()");
    DfGridFreeXC::buildFxc_LDA_method<DfOverlapX_Parallel,
                                      DfCD_Parallel,
                                      TlSymmetricMatrix,
                                      TlMatrix>();
}

void DfGridFreeXC_Parallel::buildFxc_LDA_ScaLAPACK()
{
    this->log_.info("DfGridFreeXC_Parallel::buildFxc_LDA_ScaLAPACK()");
    DfGridFreeXC::buildFxc_LDA_method<DfOverlapX_Parallel,
                                      DfCD_Parallel,
                                      TlDistributeSymmetricMatrix,
                                      TlDistributeMatrix>();
}

void DfGridFreeXC_Parallel::buildFxc_GGA()
{
    this->log_.info("DfGridFreeXC_Parallel::buildFxc_GGA()");
#ifdef HAVE_SCALAPACK
    if (this->m_bUsingSCALAPACK == true) {
        this->buildFxc_GGA_ScaLAPACK();
    } else {
        this->buildFxc_GGA_LAPACK();
    }
#else
    {
        this->buildFxc_GGA_LAPACK();
    }
#endif // HAVE_SCALAPACK
}

void DfGridFreeXC_Parallel::buildFxc_GGA_LAPACK()
{
    this->log_.info("DfGridFreeXC_Parallel::buildFxc_GGA_LAPACK()");
    DfGridFreeXC::buildFxc_GGA_method<DfOverlapX_Parallel,
                                      DfCD_Parallel,
                                      TlSymmetricMatrix,
                                      TlMatrix>();
}

void DfGridFreeXC_Parallel::buildFxc_GGA_ScaLAPACK()
{
    this->log_.info("DfGridFreeXC_Parallel::buildFxc_GGA_ScaLAPACK()");
    DfGridFreeXC::buildFxc_GGA_method<DfOverlapX_Parallel,
                                      DfCD_Parallel,
                                      TlDistributeSymmetricMatrix,
                                      TlDistributeMatrix>();
}

