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

#ifndef TLROTATELCAO_H
#define TLROTATELCAO_H

#include "TlMatrix.h"
#include "TlOrbitalInfo.h"

class TlRotateLCAO {
 public:
  TlRotateLCAO(const TlOrbitalInfo& orbInfo);
  ~TlRotateLCAO();

 public:
  TlMatrix exec(const TlMatrix& lcao, const TlMatrix& rot);

 protected:
  void rotateLCAO_typeS(const TlMatrix& lcao, const TlMatrix& rot,
                        const int row, const int maxCol, TlMatrix* ioMatrix);
  void rotateLCAO_typeP(const TlMatrix& lcao, const TlMatrix& rot,
                        const int row, const int maxCol, TlMatrix* ioMatrix);
  void rotateLCAO_typeD(const TlMatrix& lcao, const TlMatrix& rot,
                        const int row, const int maxCol, TlMatrix* ioMatrix);

 protected:
  TlOrbitalInfo orbInfo_;
};

#endif  // TLROTATELCAO_H
