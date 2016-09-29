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
#include "TlEspField.h"
#include "DfHpqX.h"
#include "Fl_Geometry.h"
#include "TlSymmetricMatrix.h"
#include "TlPosition.h"
#include "TlUtils.h"
#include "TlSerializeData.h"

#define AU_PER_ANG 1.889762

TlEspField::TlEspField(const TlSerializeData& param)
    : param_(param) {
}


TlEspField::~TlEspField()
{
}


std::vector<double> TlEspField::makeEspFld(const TlSymmetricMatrix& P,
                                           const std::vector<TlPosition>& grids)
{
    DfHpqX dfHpq(&this->param_);

    const std::size_t numOfGrids = grids.size();
    std::vector<double> values(numOfGrids);

    // electron part
    values = dfHpq.getESP(P, grids);

    // Coulomb potential generated by nuclear
    const Fl_Geometry flGeom(this->param_["coordinates"]);
    const std::size_t numOfAtoms = flGeom.getNumOfAtoms();
    for (std::size_t atomIndex = 0; atomIndex < numOfAtoms; ++atomIndex) {
        const std::string atomSymbol = flGeom.getAtomSymbol(atomIndex);
        if (atomSymbol == "X") {
            // not calculate in case of dummy charge
            continue;
        }

        const TlPosition pos = flGeom.getCoordinate(atomIndex);
        const double charge = flGeom.getCharge(atomIndex);

#pragma omp parallel for schedule(runtime)
        for (std::size_t gridIndex = 0; gridIndex < numOfGrids; ++gridIndex) {
            const TlPosition grid = grids[gridIndex];
            const double distance = pos.distanceFrom(grid);
            const double esp = charge / distance;

#pragma omp critical(TlEspField__makeEspFld)
            {
                values[gridIndex] += esp;
            }
        }
    }

    return values;
}


// atom only
std::vector<double> TlEspField::makeEspFld(const std::vector<TlPosition>& grids)
{
    const std::size_t numOfGrids = grids.size();
    std::vector<double> values(numOfGrids);

    // Coulomb potential generated by nuclear
    const Fl_Geometry flGeom(this->param_["coordinates"]);
    const std::size_t numOfAtoms = flGeom.getNumOfAtoms();
    for (std::size_t atomIndex = 0; atomIndex < numOfAtoms; ++atomIndex) {
        const std::string atomSymbol = flGeom.getAtomSymbol(atomIndex);
        if (atomSymbol == "X") {
            // not calculate in case of dummy charge
            continue;
        }

        const TlPosition pos = flGeom.getCoordinate(atomIndex);
        const double charge = flGeom.getCharge(atomIndex);

#pragma omp parallel for schedule(runtime)
        for (std::size_t gridIndex = 0; gridIndex < numOfGrids; ++gridIndex) {
            const TlPosition grid = grids[gridIndex];
            const double distance = pos.distanceFrom(grid);
            const double esp = charge / distance;

#pragma omp critical(TlEspField__makeEspFld)
            {
                values[gridIndex] += esp;
            }
        }
    }

    return values;
}

