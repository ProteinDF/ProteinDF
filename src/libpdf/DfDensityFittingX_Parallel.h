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

#ifndef DFDENSITYFITTINGX_PARALLEL_H
#define DFDENSITYFITTINGX_PARALLEL_H

#include "DfDensityFittingObject.h"
#include "DfEriX_Parallel.h"
#include "tl_dense_vector_blacs.h"

// LAPACK版
class DfDensityFittingX_Parallel
    : public DfDensityFittingTmpl<TlDenseSymmetricMatrix_BLAS_Old, TlVector_BLAS,
                                  DfEriX_Parallel> {
 public:
  DfDensityFittingX_Parallel(TlSerializeData* pPdfParam);
  virtual ~DfDensityFittingX_Parallel();

  void exec();

 protected:
  virtual TlVector_BLAS getNalpha();
  virtual TlDenseSymmetricMatrix_BLAS_Old getSabinv();
  virtual TlVector_BLAS calcTAlpha_DIRECT(const TlDenseSymmetricMatrix_BLAS_Old& P);
  virtual TlVector_BLAS getTalpha(RUN_TYPE runType, const int iteration);
  virtual void getTalpha_ROKS(TlVector_BLAS* pT_alphaA,
                              TlVector_BLAS* pT_alphaB);
  virtual TlDenseSymmetricMatrix_BLAS_Old getDiffDensityMatrix(RUN_TYPE runType);
  virtual TlDenseSymmetricMatrix_BLAS_Old getP1pq(const int nIteration);
  virtual TlDenseSymmetricMatrix_BLAS_Old getP2pq(const int nIteration);
  virtual double getLamda(const TlVector_BLAS& SabinvN,
                          const TlVector_BLAS& t_alpha,
                          const TlVector_BLAS& N_alpha,
                          const double dNumOfElec);
  virtual void saveRho(const TlVector_BLAS& rRho, RUN_TYPE runType);
};

#endif  // DFDENSITYFITTINGX_PARALLEL_H
