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

#ifndef DFGRIDFREEXC_PARALLEL_H
#define DFGRIDFREEXC_PARALLEL_H

#include "DfGridFreeXC.h"

class DfGridFreeXC_Parallel : public DfGridFreeXC {
 public:
  DfGridFreeXC_Parallel(TlSerializeData* pPdfParam);
  virtual ~DfGridFreeXC_Parallel();

 public:
  virtual void preprocessBeforeSCF();

 protected:
  virtual DfOverlapX* getDfOverlapObject();
  virtual DfXMatrix* getDfXMatrixObject();

  void preprocessBeforeSCF_LAPACK();
  void preprocessBeforeSCF_ScaLAPACK();

 protected:
  virtual void buildFxc_LDA();
  void buildFxc_LDA_LAPACK();
  void buildFxc_LDA_ScaLAPACK();

  virtual void buildFxc_GGA();
  void buildFxc_GGA_LAPACK();
  void buildFxc_GGA_ScaLAPACK();
};

#endif  // DFGRIDFREEXC_PARALLEL_H
