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

#ifndef DFDENSITYFITTINGX_H
#define DFDENSITYFITTINGX_H

#include "DfDensityFittingObject.h"
#include "DfEriX.h"
#include "DfObject.h"
#include "tl_dense_symmetric_matrix_lapack.h"
#include "tl_dense_vector_lapack.h"

class DfDensityFittingX
    : public DfDensityFittingTmpl<TlDenseSymmetricMatrix_Lapack,
                                  TlDenseVector_Lapack, DfEriX> {
   public:
    DfDensityFittingX(TlSerializeData* pPdfParam);
    virtual ~DfDensityFittingX();

   public:
    void exec();
};

#endif  // DFDENSITYFITTINGX_H
