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

#include <iostream>
#include <cmath>

#include "DfPopulation.h"
#include "Fl_Geometry.h"
#include "TlUtils.h"
#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"

DfPopulation::DfPopulation(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam), orbitalInfo_((*pPdfParam)["coordinates"],
                                        (*pPdfParam)["basis_sets"])
{
}


DfPopulation::~DfPopulation()
{
}


void DfPopulation::sumOfElectrons(int iteration, double* pAlpha, double* pBeta)
{
    this->calcPop(iteration);
    if (pAlpha != NULL) {
        *pAlpha = this->grossAtomPopA_.sum();
    }
    if (pBeta != NULL) {
        *pBeta = this->grossAtomPopB_.sum();
    }
}


TlMatrix DfPopulation::getAtomPopData(const int iteration)
{
    this->calcPop(iteration);

    const Fl_Geometry flGeom((*(this->pPdfParam_))["coordinates"]);

    TlMatrix answer;
    switch (this->m_nMethodType) {
    case METHOD_RKS:
        {
            const std::size_t dim = this->grossAtomPopA_.getSize();
            answer.resize(dim, 1);
            for (std::size_t atomIndex = 0; atomIndex < dim; ++atomIndex) {
                const double nucCharge = flGeom.getCharge(atomIndex);
                answer.set(atomIndex, 0, nucCharge - this->grossAtomPopA_.get(atomIndex));
            }
        }
        break;

    case METHOD_UKS:
        {
            const std::size_t dim = this->grossAtomPopA_.getSize();
            assert(dim == this->grossAtomPopB_.getSize());
            answer.resize(dim, 2);
            for (std::size_t atomIndex = 0; atomIndex < dim; ++atomIndex) {
                const double nucCharge = flGeom.getCharge(atomIndex);
                answer.set(atomIndex, 0, nucCharge - this->grossAtomPopA_.get(atomIndex));
                answer.set(atomIndex, 1, nucCharge - this->grossAtomPopB_.get(atomIndex));
            }
        }
        break;

    case METHOD_ROKS:
        {
            const std::size_t dim = this->grossAtomPopA_.getSize();
            answer.resize(dim, 1);
            for (std::size_t atomIndex = 0; atomIndex < dim; ++atomIndex) {
                const double nucCharge = flGeom.getCharge(atomIndex);
                answer.set(atomIndex, 0, nucCharge - this->grossAtomPopA_.get(atomIndex));
            }
        }
        break;

    default:
        break;
    }

    return answer;
}


void DfPopulation::calcPop(const int iteration)
{
    this->calcPop<TlSymmetricMatrix>(iteration);
}


double DfPopulation::getNucleiCharge()
{
    const Fl_Geometry geom((*this->pPdfParam_)["coordinates"]);

    const int numOfAtoms = this->m_nNumOfAtoms;
    double charge = 0.0;
    for (int i = 0; i < numOfAtoms; ++i) {
        charge += geom.getCharge(i);
    }

    return charge;
}


// Mulliken Analysis (Gross Atom Population)
TlVector DfPopulation::getGrossAtomPop(const TlVector& trPS)
{
    std::string output = "";

    const index_type numOfAtoms = this->m_nNumOfAtoms;
    const index_type numOfAOs = this->m_nNumOfAOs;
    
    std::vector<double> answer(numOfAtoms, 0.0);
    
#pragma omp parallel for
    for (index_type aoIndex= 0; aoIndex < numOfAOs; ++aoIndex) {
        const index_type atomIndex = this->orbitalInfo_.getAtomIndex(aoIndex);
#pragma omp atomic
        answer[atomIndex] += trPS.get(aoIndex);
    }

    return TlVector(answer);
}
    


