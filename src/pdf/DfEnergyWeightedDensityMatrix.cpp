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

#include "DfEnergyWeightedDensityMatrix.h"

DfEnergyWeightedDensityMatrix::DfEnergyWeightedDensityMatrix()
    : DfObject() {
}


DfEnergyWeightedDensityMatrix::~DfEnergyWeightedDensityMatrix()
{
}


void DfEnergyWeightedDensityMatrix::calc(RUN_TYPE runType)
{
    const int itr = this->numOfIterations_;
    const int numOfMOs = this->getNumOfMOs_;

    TlVector eps;
    eps.load(this->getOccupationPath(runType));
    assert(occ.getSize() == numOfMOs);
    {
        TlVector eig;
        eig.load(this->getEigvalPath(runType, itr));
        assert(eig.getSize() == numOfMOs);

        eps.dot(eig);
    }

    TlMatrix C = this->getCMatrix(runType, itr);
    TlMatrix Ct = C;
    Ct.transpose();
    
    TlMatrix W = Ct * eps * C;
    return W;
}



