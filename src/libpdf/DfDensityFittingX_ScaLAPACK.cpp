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

#include "DfDensityFittingX_ScaLAPACK.h"
#include "tl_dense_symmetric_matrix_blas_old.h"

// SCALAPACK ===========================================================
DfDensityFittingX_ScaLAPACK::DfDensityFittingX_ScaLAPACK(
    TlSerializeData* pPdfParam)
    : DfDensityFittingTmpl<TlDenseSymmetricMatrix_blacs, TlDistributedVector,
                           DfEriX_Parallel>(pPdfParam) {}

DfDensityFittingX_ScaLAPACK::~DfDensityFittingX_ScaLAPACK() {}

void DfDensityFittingX_ScaLAPACK::exec() { this->calc(); }
