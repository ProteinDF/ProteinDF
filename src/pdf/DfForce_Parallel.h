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

#ifndef DFFORCE_PARALLEL_H
#define DFFORCE_PARALLEL_H

#include "DfForce.h"

class DfForce_Parallel : public DfForce {
public:
    DfForce_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfForce_Parallel();

protected:
    virtual void logger(const std::string& str) const;
    
protected:
    virtual void calcForceFromCoulomb_RIJ(const RUN_TYPE runType);
    void calcForceFromCoulomb_RIJ_DC(const RUN_TYPE runType);

    virtual void calcForceFromK(const RUN_TYPE runType);
    void calcForceFromK_DC(const RUN_TYPE runType);
};


#endif // DFFORCE_PARALLEL_H
