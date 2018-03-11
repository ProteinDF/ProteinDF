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

#ifndef DFCONVERGE_ANDERSON_PARALLEL_H
#define DFCONVERGE_ANDERSON_PARALLEL_H

#include "DfConverge_Anderson.h"

class DfConverge_Anderson_Parallel : public DfConverge_Anderson {
 public:
  DfConverge_Anderson_Parallel(TlSerializeData* pPdfParam);
  virtual ~DfConverge_Anderson_Parallel();

 protected:
  virtual void convergeRhoTilde();
  virtual void convergeKSMatrix();
  virtual void convergePMatrix();

  void convergeRhoTilde_LAPACK();
  void convergeKSMatrix_LAPACK();
  void convergePMatrix_LAPACK();

  void convergeRhoTilde_ScaLAPACK();
  void convergeKSMatrix_ScaLAPACK();
  void convergePMatrix_ScaLAPACK();
};

#endif  // DFCONVERGE_ANDERSON_PARALLEL_H
