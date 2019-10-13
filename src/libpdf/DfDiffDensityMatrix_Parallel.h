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

#ifndef DFDIFFDENSITYMATRIX_PARALLEL_H
#define DFDIFFDENSITYMATRIX_PARALLEL_H

#include "DfDiffDensityMatrix.h"

class DfDiffDensityMatrix_Parallel : public DfDiffDensityMatrix {
   public:
    DfDiffDensityMatrix_Parallel(TlSerializeData* pPdfParam);
    ~DfDiffDensityMatrix_Parallel();

   public:
    virtual void exec();

   protected:
    // void calc_usingScaLAPACK(DfObject::RUN_TYPE runType, int iteration);
};

#endif  // DFDIFFDENSITYMATRIX_PARALLEL_H
