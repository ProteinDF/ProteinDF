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

#ifndef DFFUNCTIONAL_BECKE88_H
#define DFFUNCTIONAL_BECKE88_H

#include "DfFunctional.h"
#include "tl_dense_general_matrix_lapack.h"

class DfFunctional_Becke88 : public DfFunctional_GGA {
   public:
    DfFunctional_Becke88();
    virtual ~DfFunctional_Becke88();

   public:
    virtual double getFunctional(double dRhoA, double dRhoB, double dGammaAA,
                                 double dGammaAB, double dGammaBB);
    virtual double getFunctional(double dRhoA, double dGammaAA);

    virtual void getDerivativeFunctional(
        double dRhoA, double dRhoB, double dGammaAA, double dGammaAB,
        double dGammaBB, double* pRoundF_roundRhoA, double* pRoundF_roundRhoB,
        double* pRoundF_roundGammaAA, double* pRoundF_roundGammaAB,
        double* pRoundF_roundGammaBB);
    virtual void getDerivativeFunctional(double dRhoA, double dGammaAA,
                                         double* pRoundF_roundRhoA,
                                         double* pRoundF_roundGammaAA,
                                         double* pRoundF_roundGammaAB);

   public:
    // for Grid-Free ===================================================
    virtual TlDenseGeneralMatrix_Lapack getFunctionalCore(const double rhoA,
                                                          const double rhoB,
                                                          const double xA,
                                                          const double xB);

    virtual TlDenseGeneralMatrix_Lapack getDerivativeFunctionalCore(
        const double rhoA, const double rhoB, const double xA, const double xB);

   protected:
    double x(double dRho, double dGamma);

    virtual double g(double x);
    virtual double g_prime(double x);

   protected:
    static const double EPS;
    static const double INV_3;
    static const double M_4_3;
    static const double G_PARAM;
    static const double BECKE_B;

   private:
    friend class DfFunctional_B88LYP;
};

#endif  // DFFUNCTIONAL_BECKE88_H
