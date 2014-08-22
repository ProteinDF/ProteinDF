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

#include <vector>
#include <algorithm>
#include "TlDensField.h"
#include "DfCalcGridX.h"
#include "Fl_Geometry.h"
#include "TlSymmetricMatrix.h"
#include "TlPosition.h"
#include "TlUtils.h"
#include "TlSerializeData.h"

#define AU_PER_ANG 1.889762

TlDensField::TlDensField(const TlSerializeData& param)
    : param_(param) {
}


TlDensField::~TlDensField()
{
}


std::vector<double> TlDensField::makeDensFld(const TlSymmetricMatrix& P,
                              const std::vector<TlPosition>& grids)
{
    const std::size_t numOfGrids = grids.size();
    std::vector<double> values(numOfGrids, 0.0);

    DfCalcGridX dfCalcGrid(&this->param_);

#pragma omp parallel for schedule(runtime)
    for (std::size_t gridIndex = 0; gridIndex < numOfGrids; ++gridIndex) {
        const TlPosition grid = grids[gridIndex];
        double rho = 0.0;
        dfCalcGrid.gridDensity(P, grid, &rho);

#pragma omp atomic
        values[gridIndex] += rho;
    }

    return values;
}


