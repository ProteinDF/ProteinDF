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

#ifndef DFXCFUNCTIONAL_PARALLEL_H
#define DFXCFUNCTIONAL_PARALLEL_H

#include "DfXCFunctional.h"
#include "DfCalcGridX_Parallel.h"
#include "CnError.h"

// 交換・相関項の計算を行う
class DfXCFunctional_Parallel : public DfXCFunctional {
public:
    DfXCFunctional_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfXCFunctional_Parallel();

public:
    virtual void buildXcMatrix();

    virtual double getGrimmeDispersionEnergy();
    
protected:
    virtual void logger(const std::string& str) const;

    void buildXC_LAPACK();
    void buildXC_ScaLAPACK();

protected:
    virtual DfEriX* getDfEriXObject();
};


#endif // DFXCFUNCTIONAL_PARALLEL_H
