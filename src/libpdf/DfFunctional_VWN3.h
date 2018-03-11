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

#ifndef DFFUNCTIONAL_VWN3_H
#define DFFUNCTIONAL_VWN3_H

#include "DfFunctional_VWN.h"

// VWN3
class DfFunctional_VWN3 : public DfFunctional_VWN {
 public:
  DfFunctional_VWN3();
  virtual ~DfFunctional_VWN3();

 protected:
  virtual double epsilonC_PARA(double x);
  virtual double epsilonC_FERR(double x);

  virtual double epsilonCPrime_PARA(double x);
  virtual double epsilonCPrime_FERR(double x);

 protected:
  static const double VWN3_A_PARA;
  static const double VWN3_B_PARA;
  static const double VWN3_C_PARA;
  static const double VWN3_X0_PARA;
  static const double VWN3_A_FERR;
  static const double VWN3_B_FERR;
  static const double VWN3_C_FERR;
  static const double VWN3_X0_FERR;
};

#endif  // DFFUNCTIONAL_VWN3_H
