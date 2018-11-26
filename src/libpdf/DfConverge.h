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

#ifndef DFCONVERGE_H
#define DFCONVERGE_H

#include "DfObject.h"

class DfConverge : public DfObject {
 protected:
  enum ConvergeTarget { RHO_TILDE, KS_MATRIX, DENSITY_MATRIX };

 public:
  DfConverge(TlSerializeData* pPdfParam);
  virtual ~DfConverge();

 public:
  void doConverge();

 protected:
  virtual void convergeRhoTilde() = 0;
  virtual void convergeKSMatrix() = 0;
  virtual void convergePMatrix() = 0;

 protected:
  ConvergeTarget m_nConvergeTarget;
};

#endif  // DFCONVERGE
