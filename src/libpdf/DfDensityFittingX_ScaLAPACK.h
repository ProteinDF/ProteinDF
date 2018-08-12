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

#ifndef DFDENSITYFITTINGX_SCALAPACK_H
#define DFDENSITYFITTINGX_SCALAPACK_H

#include "DfDensityFittingObject.h"
#include "DfEriX_Parallel.h"
#include "tl_dense_symmetric_matrix_scalapack.h"
#include "tl_dense_vector_scalapack.h"

// 特殊化
// template<>
// TlDistributedVector DfDensityFittingTmpl<TlDenseSymmetricMatrix_Scalapack,
// TlDistributedVector, DfEriX_Parallel>::calcTAlpha_DIRECT (const
// TlDenseSymmetricMatrix_Scalapack& P)
// {
//     TlDistributedVector t_alpha(this->m_nNumOfAux);

//     DfEriX_Parallel dfEri(this->pPdfParam_);
//     dfEri.getJ_D(P, &t_alpha);

//     return t_alpha;
// }

// ScaLAPACK版
class DfDensityFittingX_ScaLAPACK
    : public DfDensityFittingTmpl<TlDenseSymmetricMatrix_Scalapack,
                                  TlDenseVector_Scalapack, DfEriX_Parallel> {
 public:
  DfDensityFittingX_ScaLAPACK(TlSerializeData* pPdfParam);
  virtual ~DfDensityFittingX_ScaLAPACK();

  void exec();
};

#endif  // DFDENSITYFITTINGX_SCALAPACK_H
