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

#ifndef DFFUNCTIONAL_SVWN_H
#define DFFUNCTIONAL_SVWN_H

#include "DfFunctional.h"
#include "DfFunctional_Slater.h"
#include "DfFunctional_VWN.h"

class DfFunctional_SVWN : public DfFunctional_LDA {
public:
    DfFunctional_SVWN();
    virtual ~DfFunctional_SVWN();

public:
    virtual double getFunctional(double dRhoA, double dRhoB);
    virtual void getDerivativeFunctional(double dRhoA, double dRhoB,
                                         double* pRoundFunctional_roundRhoA, double* pRoundFunctional_roundRhoB);

protected:
    DfFunctional_Slater m_Slater;
    DfFunctional_VWN m_VWN;
};

#endif // DFFUNCTIONAL_SVWN_H
