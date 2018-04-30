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

#ifndef DFCLEANUP_H
#define DFCLEANUP_H

#include "DfObject.h"

class DfCleanup : public DfObject {
 public:
  DfCleanup(TlSerializeData* pPdfParam);
  virtual ~DfCleanup();

 public:
  virtual void cleanup();

 protected:
  void cleanup(RUN_TYPE runType, int iteration);
  void cleanupFxc(RUN_TYPE runType, int iteration);
  void cleanupHFx(RUN_TYPE runType, int iteration);
  void cleanupFpq(RUN_TYPE runType, int iteration);
  void cleanupFprime(RUN_TYPE runType, int iteration);
  void cleanupCprime(RUN_TYPE runType, int iteration);
  void cleanupC(RUN_TYPE runType, int iteration);
  void cleanupP(RUN_TYPE runType, int iteration);
};

#endif  // DFCLEANUP_H
