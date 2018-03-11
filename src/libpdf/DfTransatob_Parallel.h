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

#ifndef DFTRANSATOB_PARALLEL_H
#define DFTRANSATOB_PARALLEL_H

#include "DfTransatob.h"

class DfTransatob_Parallel : public DfTransatob {
 public:
  DfTransatob_Parallel(TlSerializeData* pPdfParam);
  virtual ~DfTransatob_Parallel();

  virtual void DfTrsatobMain();
  virtual void DfTrsatobQclo(const std::string& fragname, int norbcut);

 protected:
  void DfTrsatobMain_SCALAPACK();
  void DfTrsatobQclo_SCALAPACK(const std::string& fragname, int norbcut);

  virtual void logger(const std::string& str) const;
};

#endif  // DFTRANSATOB_PARALLEL_H
