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

#ifndef DFTRANSFMATRIX_PARALLEL_H
#define DFTRANSFMATRIX_PARALLEL_H

#include "DfTransFmatrix.h"

class DfTransFmatrix_Parallel : public DfTransFmatrix {
 public:
  DfTransFmatrix_Parallel(TlSerializeData* pPdfParam, bool bExecDiis);
  ~DfTransFmatrix_Parallel();

 public:
  virtual void DfTrsFmatMain();
  virtual void DfTrsFmatQclo(const std::string& fragname, int norbcut);

 protected:
  virtual void logger(const std::string& str) const;

  void DfTrsFmatMain_SCALAPACK();
  void DfTrsFmatQclo_SCALAPACK(const std::string& fragname, int norbcut);
};

#endif  // DFTRANSFMATRIX_PARALLEL_H
