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

#ifndef DFPRESCF_PARALLEL_H
#define DFPRESCF_PARALLEL_H

#include "DfPreScf.h"

/// DfScfクラスの前処理を行うクラス
class DfPreScf_Parallel : public DfPreScf {
 public:
  DfPreScf_Parallel(TlSerializeData* pPdfParam);
  virtual ~DfPreScf_Parallel();

 protected:
  virtual void createInitialGuessUsingLCAO(const RUN_TYPE runType);
  virtual void createOccupation(const RUN_TYPE runType);

 protected:
  void createInitialGuessUsingLCAO_onScaLAPACK(const RUN_TYPE runType);
  void createInitialGuessUsingLCAO_onDisk(const RUN_TYPE runType);

  TlDistributeMatrix getLCAO_onScaLAPACK(const RUN_TYPE runType);

  virtual void logger(const std::string& str) const;
};

#endif  // DFPRESCF_PARALLEL_H
