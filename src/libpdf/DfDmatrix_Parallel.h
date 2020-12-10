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

#ifndef DFDMATRIX_PARALLEL_H
#define DFDMATRIX_PARALLEL_H

#include "DfDmatrix.h"
#include "tl_dense_vector_lapack.h"

class DfDmatrix_Parallel : public DfDmatrix {
   public:
    DfDmatrix_Parallel(TlSerializeData* pPdfParam);
    virtual ~DfDmatrix_Parallel();

   public:
    virtual void run();

   protected:
    void run_Scalapack();

   protected:
    virtual void checkOccupation(const TlDenseVector_Lapack& prevOcc,
                                 const TlDenseVector_Lapack& currOcc);
    virtual void printOccupation(const TlDenseVector_Lapack& occ);

   private:
    using DfDmatrix::checkOccupation;
    using DfDmatrix::printOccupation;

   protected:
    virtual TlDenseVector_Lapack getOccVtr(DfObject::RUN_TYPE runType);
};

#endif  // DFDMATRIX_PARALLEL_H
