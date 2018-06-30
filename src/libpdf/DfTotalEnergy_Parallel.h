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

#ifndef DFTOTALENERGY_PARALLEL_H
#define DFTOTALENERGY_PARALLEL_H

#include "DfTotalEnergy.h"
#include "tl_dense_vector_blacs.h"

class DfTotalEnergy_Parallel : public DfTotalEnergy {
 public:
  DfTotalEnergy_Parallel(TlSerializeData* pPdfParam);
  virtual ~DfTotalEnergy_Parallel();

 public:
  virtual void exec();

  /// ダミー原子を除いた時の全エネルギー計算を行う
  virtual void calculate_real_energy();

 protected:
  void logger(const std::string& str) const;

 protected:
  void exec_LAPACK();
  void exec_ScaLAPACK();

 protected:
  virtual void output();

 protected:
  bool m_bUseDistributeMatrix;
};

#endif  // DFTOTALENERGY_PARALLEL_H
