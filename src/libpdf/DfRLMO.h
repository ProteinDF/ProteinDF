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

#ifndef DFRLMO_H
#define DFRLMO_H

#include <vector>
#include "DfObject.h"
#include "TlMatrix.h"

class DfRLMO : public DfObject {
public:
    DfRLMO(TlSerializeData* pPdfParam);
    virtual ~DfRLMO();
    
public:
    void exec(const std::vector<index_type>& startBlockAOs);

protected:
    TlSymmetricMatrix getX();

    TlMatrix getT(const TlSymmetricMatrix& D,
                  const std::vector<index_type>& startBlockAOs);
};

#endif // DFRLMO_H
