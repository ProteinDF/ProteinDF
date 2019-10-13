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

#include <cmath>
#include <iostream>

#include "DfPopulation.h"
#include "Fl_Geometry.h"
#include "TlUtils.h"
#include "tl_dense_symmetric_matrix_lapack.h"
#include "tl_dense_vector_lapack.h"

DfPopulation::DfPopulation(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam),
      orbitalInfo_((*pPdfParam)["coordinates"], (*pPdfParam)["basis_set"]) {
    this->setNucleiCharges();
}

DfPopulation::~DfPopulation() {}

void DfPopulation::exec(const int iteration) {
    this->getAtomPopulation<TlDenseSymmetricMatrix_Lapack,
                            TlDenseVector_Lapack>(iteration);
}

// TlDenseGeneralMatrix_Lapack DfPopulation::getAtomPopData(const int iteration)
// {
//   this->calcPop(iteration);
//
//   const Fl_Geometry flGeom((*(this->pPdfParam_))["coordinates"]);
//
//   TlDenseGeneralMatrix_Lapack answer;
//   switch (this->m_nMethodType) {
//     case METHOD_RKS: {
//       const std::size_t dim = this->grossAtomPopA_.size();
//       answer.resize(dim, 1);
//       for (std::size_t atomIndex = 0; atomIndex < dim; ++atomIndex) {
//         const double nucCharge = flGeom.getCharge(atomIndex);
//         answer.set(atomIndex, 0, nucCharge -
//         this->grossAtomPopA_[atomIndex]);
//       }
//     } break;
//
//     case METHOD_UKS: {
//       const std::size_t dim = this->grossAtomPopA_.size();
//       assert(dim ==
//       static_cast<std::size_t>(this->grossAtomPopB_.getSize()));
//       answer.resize(dim, 2);
//       for (std::size_t atomIndex = 0; atomIndex < dim; ++atomIndex) {
//         const double nucCharge = flGeom.getCharge(atomIndex);
//         answer.set(atomIndex, 0, nucCharge -
//         this->grossAtomPopA_[atomIndex]); answer.set(atomIndex, 1, nucCharge
//         - this->grossAtomPopB_[atomIndex]);
//       }
//     } break;
//
//     case METHOD_ROKS: {
//       const std::size_t dim = this->grossAtomPopA_.size();
//       answer.resize(dim, 1);
//       for (std::size_t atomIndex = 0; atomIndex < dim; ++atomIndex) {
//         const double nucCharge = flGeom.getCharge(atomIndex);
//         answer.set(atomIndex, 0, nucCharge -
//         this->grossAtomPopA_[atomIndex]);
//       }
//     } break;
//
//     default:
//       break;
//   }
//
//   return answer;
// }

void DfPopulation::calcPop(const int iteration) {
    this->calcPop<TlDenseSymmetricMatrix_Lapack>(iteration);
}

void DfPopulation::setNucleiCharges() {
    const Fl_Geometry geom((*this->pPdfParam_)["coordinates"]);

    const int numOfAtoms = this->m_nNumOfAtoms;
    this->nucleiCharges_.resize(numOfAtoms);

    for (int i = 0; i < numOfAtoms; ++i) {
        this->nucleiCharges_[i] = geom.getCharge(i);
    }
}

double DfPopulation::getSumOfNucleiCharges() const {
    return this->nucleiCharges_.sum();
}

// Mulliken Analysis (Gross Atom Population)
std::valarray<double> DfPopulation::getGrossAtomPop(
    const std::valarray<double>& grossOrbPop) {
    std::string output = "";

    const index_type numOfAtoms = this->m_nNumOfAtoms;
    const index_type numOfAOs = this->m_nNumOfAOs;

    std::valarray<double> answer(0.0, numOfAtoms);

#pragma omp parallel for
    for (index_type aoIndex = 0; aoIndex < numOfAOs; ++aoIndex) {
        const index_type atomIndex = this->orbitalInfo_.getAtomIndex(aoIndex);

#pragma omp critical(DfPopulation__getGrossAtomPop)
        { answer[atomIndex] += grossOrbPop[aoIndex]; }
    }

    return answer;
}

double DfPopulation::getCharge(int atomIndex) {
    const Fl_Geometry flGeom((*(this->pPdfParam_))["coordinates"]);
    const double nucCharge = flGeom.getCharge(atomIndex);
    const double grossAtomPop = this->grossAtomPopA_[atomIndex];
    return nucCharge - grossAtomPop;
}
