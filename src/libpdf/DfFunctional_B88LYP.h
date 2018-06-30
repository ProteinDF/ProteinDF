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

#ifndef DFFUNCTIONAL_B88LYP_H
#define DFFUNCTIONAL_B88LYP_H

#include "DfFunctional.h"
#include "DfFunctional_Becke88.h"
#include "DfFunctional_LYP.h"

class DfFunctional_B88LYP : public DfFunctional_GGA {
 public:
  DfFunctional_B88LYP();
  virtual ~DfFunctional_B88LYP();

 public:
  virtual double getFunctional(double dRhoA, double dRhoB, double dGammaAA,
                               double dGammaAB, double dGammaBB);
  virtual void getDerivativeFunctional(
      double dRhoA, double dRhoB, double dGammaAA, double dGammaAB,
      double dGammaBB, double* pRoundF_roundRhoA, double* pRoundF_roundRhoB,
      double* pRoundF_roundGammaAA, double* pRoundF_roundGammaAB,
      double* pRoundF_roundGammaBB);
  // public:
  //   virtual double getEnergy(double dRhoA, double dRhoB, double dGammaAA,
  //   double dGammaAB, double dGammaBB);

  //   virtual void buildFock(double dRhoA, double dRhoB,
  //           double dGradRhoAx, double dGradRhoAy, double dGradRhoAz,
  //           double dGradRhoBx, double dGradRhoBy, double dGradRhoBz,
  //           const std::vector<DfCalcGridX::WFGrid>& aPhi,
  //           const std::vector<DfCalcGridX::WFGrid>& aGradPhiX,
  //           const std::vector<DfCalcGridX::WFGrid>& aGradPhiY,
  //           const std::vector<DfCalcGridX::WFGrid>& aGradPhiZ,
  //           double dWeight, TlMatrix_Symmetric& F);

  // for Grid-Free =======================================================
 protected:
  virtual TlVector_BLAS getFunctionalTermCoef_GF();
  virtual TlVector_BLAS getDerivativeFunctionalTermCoef_GF();

  virtual TlDenseGeneralMatrix_BLAS_old getFunctionalCore(const double rhoA,
                                                      const double rhoB,
                                                      const double xA,
                                                      const double xB);

  virtual TlDenseGeneralMatrix_BLAS_old getDerivativeFunctionalCore(
      const double rhoA, const double rhoB, const double xA, const double xB);

 protected:
  DfFunctional_Becke88 m_Becke88;
  DfFunctional_LYP m_LYP;
};

#endif  // DFFUNCTIONAL_B88LYP_H
