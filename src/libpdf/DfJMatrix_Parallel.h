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

#ifndef DFJMATRIX_PARALLEL_H
#define DFJMATRIX_PARALLEL_H

#include "DfJMatrix.h"

class DfJMatrix_Parallel : public DfJMatrix {
 public:
  DfJMatrix_Parallel(TlSerializeData* pPdfParam);
  virtual ~DfJMatrix_Parallel();

 protected:
  virtual void getJ_RI();
  virtual void getJ_CD();
  virtual void getJ_conventional();

  virtual void getJ_RI_local(TlSymmetricMatrix* pJ);
  void getJ_RI_distributed(TlDistributeSymmetricMatrix* pJ);

  virtual void getJ_CD_local(TlSymmetricMatrix* pJ);
  void getJ_CD_distributed(TlDistributeSymmetricMatrix* pJ);

  virtual void getJ_conventional_local(TlSymmetricMatrix* pJ);
  void getJ_conventional_distributed(TlDistributeSymmetricMatrix* pJ);

 protected:
  virtual void saveJMatrix(const TlSymmetricMatrix& J);
  virtual TlSymmetricMatrix getJMatrix(const int iteration);

  virtual TlVector getRho(const RUN_TYPE runType, const int iteration);
};

#endif  // DFJMATRIX_PARALLEL_H
