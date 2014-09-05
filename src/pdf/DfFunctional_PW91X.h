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

#ifndef DFFUNCTIONAL_PW91X_H
#define DFFUNCTIONAL_PW91X_H

#include "DfFunctional.h"

class DfFunctional_PW91X : public DfFunctional_GGA {
public:
    DfFunctional_PW91X();
    virtual ~DfFunctional_PW91X();

public:
    double getFunctional(double dRhoA, double dGammaAA);
    virtual double getFunctional(double dRhoA, double dRhoB, double dGammaAA, double dGammaAB, double dGammaBB);

    void getDerivativeFunctional(double dRhoA, double dGammaAA,
                                 double* pRoundF_roundRhoA, double* pRoundF_roundGammaAA, double* pRoundF_roundGammaAB);
    virtual void getDerivativeFunctional(double dRhoA, double dRhoB,
                                         double dGammaAA, double dGammaAB, double dGammaBB,
                                         double* pRoundF_roundRhoA, double* pRoundF_roundRhoB,
                                         double* pRoundF_roundGammaAA, double* pRoundF_roundGammaAB, double* pRoundF_roundGammaBB);

protected:
    double pw91x(double r, double s);
    double pw91x_roundRho(double r, double s);
    double pw91x_roundGamma(double r, double s);
};

#endif // DFFUNCTIONAL_PW91X_H
