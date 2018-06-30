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
#include "tl_dense_vector_blas.h"

class DfDmatrix_Parallel : public DfDmatrix {
 public:
  DfDmatrix_Parallel(TlSerializeData* pPdfParam);
  virtual ~DfDmatrix_Parallel();

 protected:
  virtual void main(DfObject::RUN_TYPE runType);
  void main_SCALAPACK(DfObject::RUN_TYPE runType);

  virtual void checkOccupation(const TlVector_BLAS& prevOcc,
                               const TlVector_BLAS& currOcc);
  virtual void printOccupation(const TlVector_BLAS& occ);

  virtual TlVector_BLAS getOccVtr(DfObject::RUN_TYPE runType);
};

#endif  // DFDMATRIX_PARALLEL_H
