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

class DfKMatrix_Parallel : public DfKMatrix {
 public:
  DfKMatrix_Parallel(TlSerializeData* pPdfParam);
  virtual ~DfKMatrix_Parallel();

 protected:
  // virtual void getK_RI();
  virtual void getK_CD();
  virtual void getK_conventional();

  // virtual void getK_RI_local(const RUN_TYPE runType, TlSymmetricMatrix* pJ);
  // void getK_RI_distributed(const RUN_TYPE runType,
  // TlDistributeSymmetricMatrix *pJ);

  virtual void getK_CD_local(const RUN_TYPE runType, TlSymmetricMatrix* pJ);
  void getK_CD_distributed(const RUN_TYPE runType,
                           TlDistributeSymmetricMatrix* pJ);

  virtual void getK_conventional_local(const RUN_TYPE runType,
                                       TlSymmetricMatrix* pJ);
  void getK_conventional_distributed(const RUN_TYPE runType,
                                     TlDistributeSymmetricMatrix* pJ);

 protected:
  // virtual TlSymmetricMatrix getPMatrix(const RUN_TYPE runType,
  //                                      const int iteration);
  // virtual TlSymmetricMatrix getDiffDensityMatrix(const RUN_TYPE runType);

  virtual TlSymmetricMatrix getKMatrix(const RUN_TYPE runType,
                                       const int iteration);
  virtual void saveKMatrix(const RUN_TYPE runType, const TlSymmetricMatrix& K);
};

#endif  // DFKMATRIX_PARALLEL_H
