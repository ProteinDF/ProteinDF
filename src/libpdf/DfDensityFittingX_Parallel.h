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
#include "tl_dense_symmetric_matrix_lapack.h"
#include "tl_dense_vector_lapack.h"

// LAPACK版
class DfDensityFittingX_Parallel
    : public DfDensityFittingTmpl<TlDenseSymmetricMatrix_Lapack,
                                  TlDenseVector_Lapack, DfEriX_Parallel> {
   public:
    DfDensityFittingX_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfDensityFittingX_Parallel();

    void exec();

   protected:
    virtual TlDenseVector_Lapack getNalpha();
    virtual TlDenseSymmetricMatrix_Lapack getSabinv();
    virtual TlDenseVector_Lapack calcTAlpha_DIRECT(
        const TlDenseSymmetricMatrix_Lapack& P);
    virtual TlDenseVector_Lapack getTalpha(RUN_TYPE runType,
                                           const int iteration);
    virtual void getTalpha_ROKS(TlDenseVector_Lapack* pT_alphaA,
                                TlDenseVector_Lapack* pT_alphaB);
    virtual TlDenseSymmetricMatrix_Lapack getDiffDensityMatrix(
        RUN_TYPE runType);
    virtual TlDenseSymmetricMatrix_Lapack getP1pq(const int nIteration);
    virtual TlDenseSymmetricMatrix_Lapack getP2pq(const int nIteration);
    virtual double getLamda(const TlDenseVector_Lapack& SabinvN,
                            const TlDenseVector_Lapack& t_alpha,
                            const TlDenseVector_Lapack& N_alpha,
                            const double dNumOfElec);
    virtual void saveRho(const TlDenseVector_Lapack& rRho, RUN_TYPE runType);
};

#endif  // DFDENSITYFITTINGX_PARALLEL_H
