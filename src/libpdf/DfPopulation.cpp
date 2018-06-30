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
#include "tl_dense_general_matrix_blas_old.h"
#include "tl_dense_symmetric_matrix_blas_old.h"

DfPopulation::DfPopulation(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam),
      orbitalInfo_((*pPdfParam)["coordinates"], (*pPdfParam)["basis_set"]) {}

DfPopulation::~DfPopulation() {}

double DfPopulation::getSumOfElectrons(const TlDenseSymmetricMatrix_BLAS_Old& P) {
  const TlVector_BLAS trPS = this->getPS(P);
  return trPS.sum();
}

TlDenseGeneralMatrix_BLAS_old DfPopulation::getAtomPopData(const int iteration) {
  this->calcPop(iteration);

  const Fl_Geometry flGeom((*(this->pPdfParam_))["coordinates"]);

  TlDenseGeneralMatrix_BLAS_old answer;
  switch (this->m_nMethodType) {
    case METHOD_RKS: {
      const std::size_t dim = this->grossAtomPopA_.getSize();
      answer.resize(dim, 1);
      for (std::size_t atomIndex = 0; atomIndex < dim; ++atomIndex) {
        const double nucCharge = flGeom.getCharge(atomIndex);
        answer.set(atomIndex, 0,
                   nucCharge - this->grossAtomPopA_.get(atomIndex));
      }
    } break;

    case METHOD_UKS: {
      const std::size_t dim = this->grossAtomPopA_.getSize();
      assert(dim == static_cast<std::size_t>(this->grossAtomPopB_.getSize()));
      answer.resize(dim, 2);
      for (std::size_t atomIndex = 0; atomIndex < dim; ++atomIndex) {
        const double nucCharge = flGeom.getCharge(atomIndex);
        answer.set(atomIndex, 0,
                   nucCharge - this->grossAtomPopA_.get(atomIndex));
        answer.set(atomIndex, 1,
                   nucCharge - this->grossAtomPopB_.get(atomIndex));
      }
    } break;

    case METHOD_ROKS: {
      const std::size_t dim = this->grossAtomPopA_.getSize();
      answer.resize(dim, 1);
      for (std::size_t atomIndex = 0; atomIndex < dim; ++atomIndex) {
        const double nucCharge = flGeom.getCharge(atomIndex);
        answer.set(atomIndex, 0,
                   nucCharge - this->grossAtomPopA_.get(atomIndex));
      }
    } break;

    default:
      break;
  }

  return answer;
}

void DfPopulation::calcPop(const int iteration) {
  this->calcPop<TlDenseSymmetricMatrix_BLAS_Old>(iteration);
}

double DfPopulation::getNucleiCharge() {
  const Fl_Geometry geom((*this->pPdfParam_)["coordinates"]);

  const int numOfAtoms = this->m_nNumOfAtoms;
  double charge = 0.0;
  for (int i = 0; i < numOfAtoms; ++i) {
    charge += geom.getCharge(i);
  }

  return charge;
}

// Mulliken Analysis (Gross Atom Population)
TlVector_BLAS DfPopulation::getGrossAtomPop(const TlVector_BLAS& trPS) {
  std::string output = "";

  const index_type numOfAtoms = this->m_nNumOfAtoms;
  const index_type numOfAOs = this->m_nNumOfAOs;

  std::vector<double> answer(numOfAtoms, 0.0);

#pragma omp parallel for
  for (index_type aoIndex = 0; aoIndex < numOfAOs; ++aoIndex) {
    const index_type atomIndex = this->orbitalInfo_.getAtomIndex(aoIndex);

#pragma omp critical(DfPopulation__getGrossAtomPop)
    { answer[atomIndex] += trPS.get(aoIndex); }
  }

  return TlVector_BLAS(answer);
}
