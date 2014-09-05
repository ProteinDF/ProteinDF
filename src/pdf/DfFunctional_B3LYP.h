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

#ifndef DFFUNCTIONAL_B3LYP_H
#define DFFUNCTIONAL_B3LYP_H

#include "DfFunctional.h"
#include "DfFunctional_Slater.h"
#include "DfFunctional_Becke88_ExceptLDA.h"
#include "DfFunctional_VWN3.h"
#include "DfFunctional_LYP.h"

class DfFunctional_B3LYP : public DfFunctional_GGA {
public:
    DfFunctional_B3LYP();
    virtual ~DfFunctional_B3LYP();

public:
    virtual double getFunctional(double dRhoA, double dRhoB, double dGammaAA, double dGammaAB, double dGammaBB);
    virtual double getFunctional(double dRhoA, double dGammaAA);

    virtual void getDerivativeFunctional(double dRhoA, double dRhoB, double dGammaAA, double dGammaAB, double dGammaBB,
                                         double* pRoundF_roundRhoA, double* pRoundF_roundRhoB,
                                         double* pRoundF_roundGammaAA, double* pRoundF_roundGammaAB, double* pRoundF_roundGammaBB);
    virtual void getDerivativeFunctional(double dRhoA, double dGammaAA,
                                         double* pRoundF_roundRhoA, double* pRoundF_roundGammaAA, double* pRoundF_roundGammaAB);

protected:
    virtual TlMatrix getFunctionalCore(const double rhoA, const double rhoB,
                                       const double xA, const double xB);
    virtual TlMatrix getDerivativeFunctionalCore(const double rhoA,
                                                 const double rhoB,
                                                 const double xA,
                                                 const double xB);

protected:
    DfFunctional_Slater m_LDA;
    DfFunctional_Becke88_ExceptLDA m_B88;
    DfFunctional_VWN3 m_VWN;
    DfFunctional_LYP m_LYP;

private:
    static const double LDA_COEF;
    static const double B88_COEF;
    static const double VWN_COEF;
    static const double LYP_COEF;
};

#endif // DFFUNCTIONAL_B3LYP_H
