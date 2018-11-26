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

#ifndef DFFUNCTIONAL_SLATER_H
#define DFFUNCTIONAL_SLATER_H

#include "DfFunctional.h"

class DfFunctional_Slater : public DfFunctional_LDA {
 public:
  DfFunctional_Slater();
  virtual ~DfFunctional_Slater();

 public:
  virtual double getFunctional(double dRhoA, double dRhoB);
  virtual double getFunctional(double dRhoA);

  virtual void getDerivativeFunctional(double dRhoA, double dRhoB,
                                       double* pRoundF_roundRhoA,
                                       double* pRoundF_roundRhoB);
  virtual void getDerivativeFunctional(double dRhoA, double* pRoundF_roundRhoA);

 public:
  virtual TlDenseGeneralMatrix_Lapack getFunctionalCore(const double rhoA,
                                                        const double rhoB);
  virtual TlDenseGeneralMatrix_Lapack getDerivativeFunctionalCore(
      const double rhoA, const double rhoB);

 private:
  static const double INV_3;
  static const double M_4_3;

  // static const double SLATER_ALPHA;
  static const double SLATER_COEF;
  static const double ROUND_SLATER_ROUND_RHO_COEF;
};

#endif  // DFFUNCTIONAL_SLATER_H
