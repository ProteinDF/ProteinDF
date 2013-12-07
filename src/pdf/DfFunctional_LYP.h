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

#ifndef DFFUNCTIONAL_LYP_H
#define DFFUNCTIONAL_LYP_H

#include "DfFunctional.h"

class DfFunctional_LYP : public DfFunctional_GGA {
public:
    DfFunctional_LYP();
    virtual ~DfFunctional_LYP();

public:
    virtual double getFunctional(double dRhoA, double dRhoB, double dGammaAA, double dGammaAB, double dGammaBB);
    virtual double getFunctional(double dRhoA, double dGammaAA);

    virtual void getDerivativeFunctional(double dRhoA, double dRhoB, double dGammaAA, double dGammaAB, double dGammaBB,
                                         double* pRoundF_roundRhoA, double* pRoundF_roundRhoB,
                                         double* pRoundF_roundGammaAA, double* pRoundF_roundGammaAB, double* pRoundF_roundGammaBB);
    virtual void getDerivativeFunctional(double dRhoA, double dGammaAA,
                                         double* pRoundF_roundRhoA, double* pRoundF_roundGammaAA, double* pRoundF_roundGammaAB);

private:
    double LYP(double dRhoA, double dRhoB, double dGammaAA, double dGammaAB, double dGammaBB);
    double LYP(double dRhoA, double dGammaAA);

    double omega(const double dRho);
    double delta(const double dRho);

    void roundLYP_roundGamma(double dRho, double dRhoA, double dRhoB,
                             double dOmega, double dDelta,
                             double& dRoundLYP_roundGammaAA, double& dRoundLYP_roundGammaAB, double& dRoundLYP_roundGammaBB);
    void roundLYP_roundGamma(double dRho, double dRhoA,
                             double dOmega, double dDelta,
                             double& dRoundLYP_roundGammaAA, double& dRoundLYP_roundGammaAB);

    // roundLYP_roundRhoAとroundLYP_roundRhoBを同時に求める
    void roundLYP_roundRhoAB(const double dRhoA, const double dRhoB,
                             double& roundLYP_roundRhoA, double& roundLYP_roundRhoB);

    void roundRoundLYP_roundRhoRoundGamma(double dRho, double dRhoA, double dRhoB,
                                          double dInvRho, double dRhoAB,
                                          double dOmega, double dOmegaPrime,
                                          double dDelta, double dDeltaPrime,
                                          double* pRoundRoundLYP_roundRhoARoundGammaAA, double* pRoundRoundLYP_roundRhoBRoundGammaAA,
                                          double* pRoundRoundLYP_roundRhoARoundGammaAB, double* pRoundRoundLYP_roundRhoBRoundGammaAB,
                                          double* pRoundRoundLYP_roundRhoARoundGammaBB, double* pRoundRoundLYP_roundRhoBRoundGammaBB);
    void roundRoundLYP_roundRhoRoundGamma(double dRho, double dRhoA,
                                          double dInvRho, double dRhoAB,
                                          double dOmega, double dOmegaPrime,
                                          double dDelta, double dDeltaPrime,
                                          double* pRoundRoundLYP_roundRhoARoundGammaAA,
                                          double* pRoundRoundLYP_roundRhoARoundGammaAB);

    double omega_prime(const double dRho, const double dOmega);
    double delta_prime(const double dRho, const double dDelta);

public:
    virtual TlMatrix getFunctionalCore(const double rhoA, const double rhoB,
                                       const double xA, const double xB);
    virtual TlMatrix getDerivativeFunctionalCore(const double rhoA,
                                                 const double rhoB,
                                                 const double xA,
                                                 const double xB);


private:
    static const double TOLERANCE;

    static const double INV_3;
    static const double INV_9;
    static const double M_5_3;
    static const double M_4_3;
    static const double M_8_3;
    static const double M_7_9;
    static const double M_11_3;
    static const double LYP_COEF;

    static const double LYP_PARAM_A;
    static const double LYP_PARAM_B;
    static const double LYP_PARAM_C;
    static const double LYP_PARAM_D;
    static const double LYP_PARAM_AB;
    static const double LYP_PARAM_DD;
};

#endif // DFFUNCTIONAL_LYP_H
