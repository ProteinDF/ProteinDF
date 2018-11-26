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

#include "TlDensField.h"
#include <algorithm>
#include <vector>
#include "DfCalcGridX.h"
#include "Fl_Geometry.h"
#include "TlPosition.h"
#include "TlSerializeData.h"
#include "TlUtils.h"
#include "tl_dense_symmetric_matrix_lapack.h"

#define AU_PER_ANG 1.889762

TlDensField::TlDensField(const TlSerializeData& param) : param_(param) {}

TlDensField::~TlDensField() {}

std::vector<double> TlDensField::makeDensFld(
    const TlDenseSymmetricMatrix_Lapack& P,
    const std::vector<TlPosition>& grids) {
  const std::size_t numOfGrids = grids.size();
  std::vector<double> values(numOfGrids, 0.0);

  DfCalcGridX dfCalcGrid(&this->param_);

#pragma omp parallel for schedule(runtime)
  for (std::size_t gridIndex = 0; gridIndex < numOfGrids; ++gridIndex) {
    const TlPosition grid = grids[gridIndex];
    double rho = 0.0;
    dfCalcGrid.gridDensity(P, grid, &rho);

#pragma omp critical(TlDensField__makeDensFld)
    { values[gridIndex] += rho; }
  }

  return values;
}
