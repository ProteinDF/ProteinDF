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
#include "TlDistributeSymmetricMatrix.h"
#include "TlDistributeVector.h"

// LAPACK版
class DfDensityFittingX_Parallel
    : public DfDensityFittingTmpl<TlSymmetricMatrix, TlVector,
                                  DfEriX_Parallel> {
 public:
  DfDensityFittingX_Parallel(TlSerializeData* pPdfParam);
  virtual ~DfDensityFittingX_Parallel();

  void exec();

 protected:
  virtual TlVector getNalpha();
  virtual TlSymmetricMatrix getSabinv();
  virtual TlVector calcTAlpha_DIRECT(const TlSymmetricMatrix& P);
  virtual TlVector getTalpha(RUN_TYPE runType, const int iteration);
  virtual void getTalpha_ROKS(TlVector* pT_alphaA, TlVector* pT_alphaB);
  virtual TlSymmetricMatrix getDiffDensityMatrix(RUN_TYPE runType);
  virtual TlSymmetricMatrix getP1pq(const int nIteration);
  virtual TlSymmetricMatrix getP2pq(const int nIteration);
  virtual double getLamda(const TlVector& SabinvN, const TlVector& t_alpha,
                          const TlVector& N_alpha, const double dNumOfElec);
  virtual void saveRho(const TlVector& rRho, RUN_TYPE runType);
};

#endif  // DFDENSITYFITTINGX_PARALLEL_H
