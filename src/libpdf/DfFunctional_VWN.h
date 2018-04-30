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

#ifndef DFFUNCTIONAL_VWN_H
#define DFFUNCTIONAL_VWN_H

#include "DfFunctional.h"

// VWN means VWN5
// if you would like to use VWN3, please call "DfFunctional_VWN3".
class DfFunctional_VWN : public DfFunctional_LDA {
 public:
  DfFunctional_VWN();
  virtual ~DfFunctional_VWN();

 public:
  virtual double getFunctional(double dRhoA, double dRhoB);
  virtual double getFunctional(double dRhoA);

  virtual void getDerivativeFunctional(double dRhoA, double dRhoB,
                                       double* pRoundF_roundRhoA,
                                       double* pRoundF_roundRhoB);
  virtual void getDerivativeFunctional(double dRhoA, double* pRoundF_roundRhoA);

 public:
  virtual TlMatrix getFunctionalCore(const double rhoA, const double rhoB);
  virtual TlMatrix getDerivativeFunctionalCore(const double rhoA,
                                               const double rhoB);

 protected:
  // for UKS
  double VWN(double dRhoA, double dRhoB);

  // for RKS
  double VWN(double dRhoA);

  // get epsilon C
  // for UKS
  double epsilonC(const double x, const double zeta);

  // get epsilon C
  // for RKS
  double epsilonC(const double x);

  double g(const double zeta);

  double epsilonC(const double A, const double b, const double c,
                  const double x0, const double x);

  double h(const double dEc_p, const double dEc_f, const double dEc_a);

  double epsilonCPrime(const double A, const double b, const double c,
                       const double x0, const double x);

  double h_prime(const double dEc_p, const double dEc_f, const double dEc_a,
                 const double dEc_p_prime, const double dEc_f_prime,
                 const double dEc_a_prime);

  double g_prime(const double zeta);

  void roundVWN_roundRho(const double dRhoA, const double dRhoB,
                         double* pRoundF_roundRhoA, double* pRoundF_roundRhoB);
  void roundVWN_roundRho(const double dRhoA, double* pRoundF_roundRhoA);

 protected:
  virtual double epsilonC_PARA(double x);
  virtual double epsilonC_FERR(double x);
  virtual double epsilonC_ANTI(double x);
  virtual double epsilonCPrime_PARA(double x);
  virtual double epsilonCPrime_FERR(double x);
  virtual double epsilonCPrime_ANTI(double x);

 protected:
  static const double TOLERANCE;

  static const double M_3_4PI;
  static const double M_3_2;
  static const double INV_3;
  static const double INV_6;
  static const double EC_COEF;
  static const double M_4_3;
  static const double M_9_8;

  static const double VWN5_A_PARA;
  static const double VWN5_B_PARA;
  static const double VWN5_C_PARA;
  static const double VWN5_X0_PARA;
  static const double VWN5_A_FERR;
  static const double VWN5_B_FERR;
  static const double VWN5_C_FERR;
  static const double VWN5_X0_FERR;
  static const double VWN_A_ANTI;
  static const double VWN_B_ANTI;
  static const double VWN_C_ANTI;
  static const double VWN_X0_ANTI;
};

#endif  // DFFUNCTIONAL_VWN_H
