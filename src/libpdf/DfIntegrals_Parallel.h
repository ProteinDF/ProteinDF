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

#ifndef DFINTEGRALS_PARALLEL_H
#define DFINTEGRALS_PARALLEL_H

#include "DfIntegrals.h"

class DfHpq;
class DfOverlap;

class DfIntegrals_Parallel : public DfIntegrals {
 public:
  explicit DfIntegrals_Parallel(TlSerializeData* param = NULL);
  virtual ~DfIntegrals_Parallel();

 protected:
  virtual void logger(const std::string& str) const;

 protected:
  virtual DfCD* getDfCDObject();
  virtual DfXMatrix* getDfXMatrixObject();
  virtual DfInvMatrix* getDfInvMatrixObject();
  virtual DfGenerateGrid* getDfGenerateGridObject();

 protected:
  virtual void createHpqMatrix();
  void createHpqMatrix_LAPACK();
  void createHpqMatrix_ScaLAPACK();

  virtual void createOverlapMatrix();
  void createOverlapMatrix_LAPACK();
  void createOverlapMatrix_ScaLAPACK();

  virtual void createERIMatrix();
  void createERIMatrix_LAPACK();
  void createERIMatrix_ScaLAPACK();

 protected:
  virtual void outputStartTitle(const std::string& stepName,
                                const char lineChar = '-');
  virtual void outputEndTitle(const std::string& stepName = "",
                              const char lineChar = '-');
};

#endif  // DFINTEGRALS_PARALLEL_H
