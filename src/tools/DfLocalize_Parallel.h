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

#ifndef DFLOCALIZE_PARALLEL_H
#define DFLOCALIZE_PARALLEL_H

#include <vector>
#include "DfLocalize.h"

class DfLocalize_Parallel : public DfLocalize {
public:
    DfLocalize_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfLocalize_Parallel();

    virtual void localize(const std::string& inputCMatrixPath = "");

protected:
    virtual int getJobItem(DfLocalize::JobItem* pJob, bool isInitialized = false);

protected:
    std::vector<bool> jobFinishedList_;
    std::vector<bool> jobOccupiedOrb_;
};

#endif // DFLOCALIZE_PARALLEL_H
