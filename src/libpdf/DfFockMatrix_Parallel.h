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

#ifndef DFFOCKMATRIX_PARALLEL_H
#define DFFOCKMATRIX_PARALLEL_H

#include "DfFockMatrix.h"
#include "tl_dense_symmetric_matrix_lapack.h"

class TlDenseSymmetricMatrix_Scalapack;

class DfFockMatrix_Parallel : public DfFockMatrix {
   public:
    DfFockMatrix_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfFockMatrix_Parallel();

   protected:
    virtual void logger(const std::string& str) const;

   protected:
    virtual void mainDIRECT_RKS();
    // void mainDIRECT_RKS_LAPACK();
    // void mainDIRECT_RKS_ScaLAPACK();

    virtual void mainDIRECT_UKS();
    virtual void mainDIRECT_ROKS();

    virtual void setXC_RI(RUN_TYPE nRunType, TlDenseSymmetricMatrix_Lapack& F);
    virtual void setXC_DIRECT(RUN_TYPE nRunType,
                              TlDenseSymmetricMatrix_Lapack& F);

    virtual void setCoulomb(METHOD_TYPE nRunType,
                            TlDenseSymmetricMatrix_Lapack& F);
    void setCoulomb(const METHOD_TYPE nMethodType,
                    TlDenseSymmetricMatrix_Scalapack& F);

    virtual TlDenseSymmetricMatrix_Lapack getFpqMatrix(RUN_TYPE nRunType,
                                                       int nIteration);

    virtual TlDenseVector_Lapack getRho(RUN_TYPE nRunType, int nIteration);
    virtual TlDenseVector_Lapack getMyu(RUN_TYPE nRunType, int nIteration);
    // virtual void saveFpqMatrix(RUN_TYPE nRunType, const
    // TlDenseSymmetricMatrix_Lapack&
    // F);
};

#endif  // DFFOCKMATRIX_PARALLEL_H
