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

#include "DfFunctional_Becke88_ExceptLDA.h"
#include "TlMath.h"

DfFunctional_Becke88_ExceptLDA::DfFunctional_Becke88_ExceptLDA() {}

DfFunctional_Becke88_ExceptLDA::~DfFunctional_Becke88_ExceptLDA() {}

double DfFunctional_Becke88_ExceptLDA::g(const double x) {
  const double bx = BECKE_B * x;
  const double arcsinhx = TlMath::arcsinh(x);

  const double dAnswer = -((bx * x) / (1.0 + 6.0 * bx * arcsinhx));

  return dAnswer;
}
