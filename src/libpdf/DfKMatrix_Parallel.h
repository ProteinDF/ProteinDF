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

#ifndef DFKMATRIX_PARALLEL_H
#define DFKMATRIX_PARALLEL_H

#include "DfKMatrix.h"
#include "tl_dense_symmetric_matrix_lapack.h"

class TlDenseSymmetricMatrix_Scalapack;

class DfKMatrix_Parallel : public DfKMatrix {
 public:
  DfKMatrix_Parallel(TlSerializeData* pPdfParam);
  virtual ~DfKMatrix_Parallel();

 protected:
  // virtual void getK_RI();
  virtual void getK_CD();
  virtual void getK_conventional();

  // virtual void getK_RI_local(const RUN_TYPE runType,
  // TlDenseSymmetricMatrix_Lapack*
  // pJ);
  // void getK_RI_distributed(const RUN_TYPE runType,
  // TlDenseSymmetricMatrix_Scalapack *pJ);

  virtual void getK_CD_replica(const RUN_TYPE runType);
  void getK_CD_distributed(const RUN_TYPE runType);

  virtual void getK_conventional_local(const RUN_TYPE runType,
                                       TlDenseSymmetricMatrix_Lapack* pJ);
  void getK_conventional_distributed(const RUN_TYPE runType,
                                     TlDenseSymmetricMatrix_Scalapack* pJ);

 protected:
  // virtual TlDenseSymmetricMatrix_Lapack getPMatrix(const RUN_TYPE runType,
  //                                      const int iteration);
  // virtual TlDenseSymmetricMatrix_Lapack getDiffDensityMatrix(const RUN_TYPE
  // runType);

  virtual TlDenseSymmetricMatrix_Lapack getKMatrix(const RUN_TYPE runType,
                                                   const int iteration);
  virtual void saveKMatrix(const RUN_TYPE runType,
                           const TlDenseSymmetricMatrix_Lapack& K);
};

#endif  // DFKMATRIX_PARALLEL_H
